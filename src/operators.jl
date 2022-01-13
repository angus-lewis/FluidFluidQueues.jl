abstract type FluidFluidOperator{N} <: AbstractSparseArray{Float64,Int,N} end

# show(io::IO, mime::MIME"text/plain", B::FluidFluidOperator) = show(io, mime, B.array)

SparseArrays.nonzeros(B::FluidFluidOperator) = nonzeros(B.array)
SparseArrays.nonzeroinds(R::FluidFluidOperator) = SparseArrays.nonzeroinds(R.array)

struct FluidFluidGenerator <: FluidFluidOperator{2}
    array::SparseMatrixCSC{Float64,Int}
    ffq::FluidFluidQueue
end
FluidFluidGenerator(ffq::FluidFluidQueue) = FluidFluidGenerator(build_full_generator(ffq.dq).B,ffq)

struct RatesOperator <: FluidFluidOperator{1} 
    array::SparseVector{Float64,Int}
    ffq::FluidFluidQueue
end

function RatesOperator(ffq::FluidFluidQueue) 
    r = abs.([ffq.rates.lwr[:]; 
        repeat(ffq.rates.interior[:],inner=n_bases_per_cell(ffq.dq));
        ffq.rates.upr[:]])
    return RatesOperator(SparseVector{Float64,Int}([r[i]!=0.0 ? 1.0./r[i] : 0.0 for i in 1:length(r)]), ffq)
end

struct InOutGenerator <: FluidFluidOperator{2}
    array::SparseMatrixCSC{Float64,Int}
    ffq::FluidFluidQueue
end

function InOutGenerator(B::FluidFluidGenerator,s::Float64) 
    R = RatesOperator(B.ffq)
    num_nz = size(B[NotZero,:],1)
    num_z = size(B,1) - num_nz
    Ds = (SparseMatrixCSC{Float64,Int}(diagm(R.array))*B.array) - s*I(num_nz+num_z) + 
        (SparseMatrixCSC{Float64,Int}(diagm(R.array))*B[:,Zero])*
        SparseMatrixCSC{Float64,Int}(inv(s*I(num_z) - Matrix(B[Zero,Zero])))*B[Zero,:]
    return InOutGenerator(Ds,B.ffq)
end

function InOutGenerator(ffq::FluidFluidQueue,s::Float64)
    B = FluidFluidGenerator(ffq)
    return InOutGenerator(B,s)
end

function build_psi(
    D::InOutGenerator; 
    MaxIters::Int = 1000, 
    err::Float64 = 1e-9,
    v::Bool = false,
)

    exitflag = ""

    Psi = zeros(Float64, size(D[Plus,Minus]))
    A = D[Plus,Plus]
    B = D[Minus,Minus]
    C = D[Plus,Minus]
    # RA, QA = LinearAlgebra.schur(Matrix(A)) # uncomment for algorithm 1
    # RB, QB = LinearAlgebra.schur(Matrix(B)) # uncomment for algorithm 1
    OldPsi = Psi
    flag = 1
    for n = 1:MaxIters
        ## line below goes with algorithm 4, is quadratically convergent
        Psi = sylvester(Matrix(A), Matrix(B), Matrix(C))
        ## line below goes with algorithm 1
        # Psi = mysyl(RA,QA,RB,QB,C)
        @show maximum(abs.(OldPsi - Psi))
        @show n
        if maximum(abs.(OldPsi - Psi)) < err
            flag = 0
            exitflag = string(
                "Reached err tolerance in ",
                n,
                " iterations with error ",
                string(maximum(abs.(OldPsi - Psi))),
            )
            break
        elseif any(isnan.(Psi))
            flag = 0
            exitflag = string("Produced NaNs at iteration ", n)
            break
        end
        OldPsi = copy(Psi)
        ## Algorithm 4 of Bean, O'Reilly, Taylor, 2008, 
        ## Algorithms for the Laplace–Stieltjes Transforms of First Return Times for Stochastic Fluid Flows, 
        ## Methodol Comput Appl Probab (2008) 10:381–408, DOI 10.1007/s11009-008-9077-3
        A .= D[Plus,Plus] + Psi * D[Minus,Plus]
        B .= D[Minus,Minus] + D[Minus,Plus] * Psi
        C .= D[Plus,Minus] - Psi * D[Minus,Plus] * Psi
        ## Algorithm 1 of the above citation
        # A, B dont change, need only comput schur decomp once 
        # C = D["+-"] + Psi * D["-+"] * Psi
    end
    if flag == 1
        exitflag = string(
            "Reached Max Iters ",
            MaxIters,
            " with error ",
            string(maximum(abs.(OldPsi - Psi))),
        )
    end
    v && println("UPDATE: Iterations for Ψ exited with flag: ", exitflag)
    return Psi
end

function build_psi(ffq::FluidFluidQueue, s::Union{Real,Complex}; kwargs...)
    return build_psi(InOutGenerator(ffq,s; kwargs...))
end

function build_xi(B::FluidFluidGenerator,Ψ::Array{Float64,2})
    # BBullet = [B["--"] B["-0"]; B["0-"] B["00"]]
    # invB = inv(Matrix(Bbullet))

    # BBulletPlus = [B["-+"]; B["0+"]]

    # solve the linear system -[ξ 0] [Bmm Bm0; B₀₋ B00]^-1 [B₋₊; B₀₊]Ψ = ξ
    # writing out the system it turns out we only need the index -- and -0
    # blocks of the inverse. Wikipedia tells us that these are

    # the next two lines are a fancy way to do an inverse of a block matrix
    tempMat = inv(Matrix(B[Zero,Zero]))
    invBmm = inv(B[Minus,Minus] - B[Minus,Zero]*tempMat*B[Zero,Minus])
    invBm0 = -invBmm*B[Minus,Zero]*tempMat

    A = -(invBmm*B[Minus,Plus]*Ψ + invBm0*B[Zero,Plus]*Ψ + I)
    b = zeros(1,size(invBmm,1))

    A[:,1] .= 1.0 # normalisation conditions
    b[1] = 1.0 # normalisation conditions

    ξ = b/A
    coeffs = zeros(size(B,1))
    coeffs[index(B.ffq,Minus)] = ξ
    return SFMDistribution(coeffs,B.ffq.dq)
end

function build_limit_dist_operators(B::FluidFluidGenerator,D::InOutGenerator,
    R::RatesOperator,Ψ::Array{<:Real,2},ξ::SFMDistribution)

    B00inv = inv(Matrix(B[Zero,Zero]))
    invBmm = inv(Matrix(D[Minus,Minus]))# inv(B[Minus,Minus] - B[Minus,Zero]*B00inv*B[Zero,Minus])
    invBm0 = -invBmm*B[Minus,Zero]*B00inv

    # B00inv = inv(Matrix(B["00"]))
    # invBmm = inv(B["--"] - B["-0"]*B00inv*B["0-"])
    # invBm0 = -invBmm*B["-0"]*B00inv

    αp = Array(transpose(ξ.coeffs[index(B.ffq,Minus)])) * -[invBmm invBm0]

    K = D[Plus,Plus] + Ψ * D[Minus,Plus]

    n₊ = size(K,1)
    n₀ = size(invBm0,2)
    n₋ = size(invBm0,1)

    # BBulletPlus = [B[("-","+"),(:,:)]; B[("0","+"),(:,:)]]

    αintegralPibullet = ((αp * [B[Minus,Plus]; B[Zero,Plus]]) / -K) * [diagm(R[Plus]) Ψ*diagm(R[Minus])]
    αintegralPi0 = -αintegralPibullet * [B[Plus,Zero];B[Minus,Zero]] * B00inv

    α = sum(αintegralPibullet) + sum(αintegralPi0) + sum(αp)

    p = αp ./ α
    integralPibullet = αintegralPibullet ./ α
    integralPi0 = αintegralPi0 ./ α

    marginalX = zeros(Float64, n₊ + n₋ + n₀)
    marginalX[index(B.ffq,Plus)] = integralPibullet[1:n₊]
    marginalX[index(B.ffq,Minus)] = integralPibullet[(n₊+1):end] + p[1:n₋]
    marginalX[index(B.ffq,Zero)] = integralPi0[:] + p[(n₋+1):end]

    Y_point_mass = zeros(Float64, n₊ + n₋ + n₀)
    Y_point_mass[index(B.ffq,Minus)] = p[1:n₋]
    Y_point_mass[index(B.ffq,Zero)] = p[n₋+1:end]
    return SFMDistribution(marginalX,B.ffq.dq), SFMDistribution(Y_point_mass,B.ffq.dq), K
end