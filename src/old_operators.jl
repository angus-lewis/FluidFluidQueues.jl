
# function MakeRDict(
#     R::Union{Array{<:Real,2},SparseArrays.SparseMatrixCSC{<:Real,Int64}},
#     model::Model, 
#     mesh::Mesh,
#     )

#     ## Make a Dictionary so that the blocks of B are easy to access
#     N₋ = sum(model.C.<=0)
#     N₊ = sum(model.C.>=0)

#     RDict = Dict{Tuple{String,Union{Int,Colon}},SparseArrays.SparseMatrixCSC{Float64,Int64}}()

#     ppositions = cumsum(model.C .<= 0)
#     qpositions = cumsum(model.C .>= 0)
#     for ℓ in ["+", "-"]
#         for i = 1:NPhases(model)
#             FilBases = repeat(mesh.Fil[(ℓ,i)]', NBases(mesh), 1)[:]
#             pitemp = falses(N₋)
#             qitemp = falses(N₊)
#             if model.C[i] <= 0
#                 pitemp[ppositions[i]] = mesh.Fil["p"*ℓ,i][1]
#             end
#             if model.C[i] >= 0
#                 qitemp[qpositions[i]] = mesh.Fil["q"*ℓ,i][1]
#             end
#             i_idx = [
#                 pitemp
#                 falses((i - 1) * TotalNBases(mesh))
#                 FilBases
#                 falses(NPhases(model) * TotalNBases(mesh) - i * TotalNBases(mesh))
#                 qitemp
#             ]
#             RDict[(ℓ,i)] = R[i_idx, i_idx]
#         end
#         FlBases =
#             [mesh.Fil["p"*ℓ,:]; repeat(mesh.Fil[ℓ,:]', NBases(mesh), 1)[:]; mesh.Fil["q"*ℓ,:]]
#         RDict[ℓ,:] = R[FlBases, FlBases]
#     end
#     return RDict
# end

"""
# Construct the DG approximation to the operator `R`.

    MakeR(
        model::Model,
        mesh::DGMesh;
        approxType::String = "projection",
    )

# Arguments
- `model`: a Model object
- `Mmesh`: a Mesh object
- `approxType::String`: (optional) either "interpolation" or
    "projection" (default).

# Output
- a tuple with keys
    - `:R::SparseArrays.SparseMatrixCSC{Float64,Int64}`: an approximation to R
        for the whole space. If ``rᵢ(x)=0`` on any cell, the corresponding
        elements of R are zero.
    - `:RDict::Dict{Tuple{String,Union{Int,Colon}},SparseArrays.SparseMatrixCSC{Float64,Int64}}`: a
        disctionary containing sub-blocks of R. Keys are of the form
        `"PhaseSign"` or just `"Sign"`. i.e. `"1-"` cells in ``Fᵢ⁻``, and
        `"-"` for cells in ``∪ᵢFᵢ⁻``.
"""
function MakeR(
    model::Model,
    mesh::Mesh;
    approxType::String = "projection",
    v::Bool = false,
)
    V = vandermonde(NBases(mesh))

    EvalR = 1.0 ./ model.r.a(CellNodes(mesh)[:])

    N₋ = sum(model.C .<= 0)
    N₊ = sum(model.C .>= 0)

    R = SparseArrays.spzeros(
        Float64,
        N₋ + N₊ + TotalNBases(mesh) * NPhases(model),
        N₋ + N₊ + TotalNBases(mesh) * NPhases(model),
    )
    # at the boundaries
    R[1:N₋, 1:N₋] = (1.0 ./ model.r.a(model.Bounds[1,1])[model.C .<= 0]).*LinearAlgebra.I(N₋)
    R[(end-N₊+1):end, (end-N₊+1):end] =  (1.0 ./ model.r.a(model.Bounds[1,end])[model.C .>= 0]).* LinearAlgebra.I(N₊)

    # on the interior
    for n = 1:(NIntervals(mesh)*NPhases(model))
        if Basis(mesh) == "legendre"
            if approxType == "interpolation"
                leftM = V.V'
                rightM = V.inv'
            elseif approxType == "projection"
                leftM = V.V' * LinearAlgebra.diagm(V.w)
                rightM = V.V
            end
            temp = leftM*LinearAlgebra.diagm(EvalR[NBases(mesh)*(n-1).+(1:NBases(mesh))])*rightM
        elseif Basis(mesh) == "lagrange"
            if approxType == "interpolation"
                temp = LinearAlgebra.diagm(EvalR[NBases(mesh)*(n-1).+(1:NBases(mesh))])
            elseif approxType == "projection"
                # the first term, LinearAlgebra.diagm(EvalR[NBases(mesh)*(n-1).+(1:NBases(mesh))])
                # is the quadrature approximation of M^r. The quadrature weights to not
                # appear since they cancel when we transform to integral/probability
                # representation. The second term V.V*V.V' is Minv. The last term
                # LinearAlgebra.diagm(V.w)is a result of the conversion to probability
                # / integral representation.
                temp = LinearAlgebra.diagm(EvalR[NBases(mesh)*(n-1).+(1:NBases(mesh))])*V.V*V.V'*LinearAlgebra.diagm(V.w)
            end
        elseif Basis(mesh) == ""
            temp = LinearAlgebra.diagm(EvalR[NBases(mesh)*(n-1).+(1:NBases(mesh))])
        end
        R[NBases(mesh)*(n-1).+(1:NBases(mesh)).+N₋, NBases(mesh)*(n-1).+(1:NBases(mesh)).+N₋] = temp
    end

    # construc the dictionary
    RDict = MakeRDict(R,model,mesh)

    out = (R=R, RDict=RDict)
    v && println("UPDATE: R object created with keys ", keys(out))
    return out
end

"""
Construct the operator `D(s)` from `B, R`.

    MakeD(
        mesh::Mesh,
        B::Generator,
        R::NamedTuple{(:R, :RDict)},
    )

# Arguments
- `R`: a tuple as constructed by MakeR
- `B`: a tuple as constructed by MakeB
- `model`: a Model object
- `mesh`: a Mesh object

# Output
- `DDict::Dict{Tuple{String,String},Function(s::Real)}`: a dictionary of functions. Keys are
  of the for `"ℓm"` where `ℓ,m∈{+,-}`. Values are functions with one argument.
  Usage is along the lines of `D["+-"](s=1)`.
"""
function MakeD(
    # mesh::Mesh,
    B::Generator,
    R::NamedTuple{(:R, :RDict)},
    v::Bool = false,
)
    DDict = Dict{Tuple{String,String},Any}()
    for ℓ in ["+", "-"], m in ["+", "-"]
        nℓ = sum(B.Fil["p"*ℓ,:]) + sum(B.Fil[ℓ,:]) * size(B.D,1) + sum(B.Fil["q"*ℓ,:])
        Idℓ = SparseArrays.sparse(LinearAlgebra.I,nℓ,nℓ)
        if any(B.Fil["p0",:]) || any(B.Fil["0",:]) || any(B.Fil["q0",:]) # in("0", model.Signs)
            n0 = sum(B.Fil["p0",:]) +
                sum(B.Fil["0",:]) * size(B.D,1) +
                sum(B.Fil["q0",:])
            Id0 = SparseArrays.sparse(LinearAlgebra.I,n0,n0)
            DDict[ℓ,m] = function (; s::Real = 0)
                return if (ℓ == m)
                    R.RDict[ℓ,:] * (
                        B[(ℓ,m),(:,:)] - s * Idℓ +
                        B[(ℓ,"0"),(:,:)] * inv(Matrix(s * Id0 - B[("0","0"),(:,:)])) * B[("0",m),(:,:)]
                    )
                else
                    R.RDict[ℓ,:] * (
                        B[(ℓ,m),(:,:)] +
                        B[(ℓ,"0"),(:,:)] * inv(Matrix(s * Id0 - B[("0","0"),(:,:)])) * B[("0",m),(:,:)]
                    )
                end
            end # end function
        else
            DDict[ℓ,m] = function (; s::Real = 0)
                return if (ℓ == m)
                    R.RDict[ℓ,:] * (B[(ℓ,m),(:,:)] - s * Idℓ)
                else
                    R.RDict[ℓ,:] * B[(ℓ,m),(:,:)]
                end
            end # end function
        end # end if ...
    end # end for ℓ ...
    v && println("UPDATE: D(s) operator created with keys ", keys(DDict))
    return (DDict = DDict)
end
# below to be deprecated
function MakeD(
    mesh::Mesh,
    B::Generator,
    R::NamedTuple{(:R, :RDict)},
    v::Bool = false,
)
    DDict = Dict{Tuple{String,String},Any}()
    for ℓ in ["+", "-"], m in ["+", "-"]
        nℓ = sum(mesh.Fil["p"*ℓ,:]) + sum(mesh.Fil[ℓ,:]) * NBases(mesh) + sum(mesh.Fil["q"*ℓ,:])
        Idℓ = SparseArrays.sparse(LinearAlgebra.I,nℓ,nℓ)
        if any(mesh.Fil["p0",:]) || any(mesh.Fil["0",:]) || any(mesh.Fil["q0",:]) # in("0", model.Signs)
            n0 = sum(mesh.Fil["p0",:]) +
                sum(mesh.Fil["0",:]) * NBases(mesh) +
                sum(mesh.Fil["q0",:])
            Id0 = SparseArrays.sparse(LinearAlgebra.I,n0,n0)
            DDict[ℓ,m] = function (; s::Real = 0)
                return if (ℓ == m)
                    R.RDict[ℓ,:] * (
                        B[(ℓ,m),(:,:)] - s * Idℓ +
                        B[(ℓ,"0"),(:,:)] * inv(Matrix(s * Id0 - B[("0","0"),(:,:)])) * B[("0",m),(:,:)]
                    )
                else
                    R.RDict[ℓ,:] * (
                        B[(ℓ,m),(:,:)] +
                        B[(ℓ,"0"),(:,:)] * inv(Matrix(s * Id0 - B[("0","0"),(:,:)])) * B[("0",m),(:,:)]
                    )
                end
            end # end function
        else
            DDict[ℓ,m] = function (; s::Real = 0)
                return if (ℓ == m)
                    R.RDict[ℓ,:] * (B[(ℓ,m),(:,:)] - s * Idℓ)
                else
                    R.RDict[ℓ,:] * B[(ℓ,m),(:,:)]
                end
            end # end function
        end # end if ...
    end # end for ℓ ...
    v && println("UPDATE: D(s) operator created with keys ", keys(DDict))
    return (DDict = DDict)
end

"""
Construct and evaluate ``Ψ(s)``.

Uses newtons method to solve the Ricatti equation
``D⁺⁻(s) + Ψ(s)D⁻⁺(s)Ψ(s) + Ψ(s)D⁻⁻(s) + D⁺⁺(s)Ψ(s) = 0.``

    PsiFun(; s = 0, D, MaxIters = 1000, err = 1e-8)

# Arguments
- `s::Real`: a value to evaluate the LST at
- `D`: a `Dict{String,Function(s::Real)}` as output from MakeD
- `MaxIters::Int`: the maximum number of iterations of newtons method
- `err::Float64`: an error tolerance for terminating newtons method. Terminates
    when `max(Ψ_{n} - Ψ{n-1}) .< eps`.

# Output
- `Ψ(s)::Array{Float64,2}`: a matrix approxiamtion to ``Ψ(s)``.
"""
function PsiFun(
    D::Dict{Tuple{String,String},Any}; 
    s::Real = 0, 
    MaxIters::Int = 1000, 
    err::Float64 = 1e-8,
    v::Bool = false,
)

    exitflag = ""

    EvalD = Dict{Tuple{String,String},SparseArrays.SparseMatrixCSC{Float64,Int64}}(("+","-") => D["+","-"](s = s))
    Dimensions = size(EvalD["+","-"])
    for ℓ in [("+","+"), ("-","-"), ("-","+")]
        EvalD[ℓ] = D[ℓ](s = s)
    end
    Psi = zeros(Float64, Dimensions)
    A = EvalD["+","+"]
    B = EvalD["-","-"]
    C = EvalD["+","-"]
    # RA, QA = LinearAlgebra.schur(Matrix(A)) # uncomment for algorithm 1
    # RB, QB = LinearAlgebra.schur(Matrix(B)) # uncomment for algorithm 1
    OldPsi = Psi
    flag = 1
    for n = 1:MaxIters
        ## line below goes with algorithm 4, is quadratically convergent
        Psi = LinearAlgebra.sylvester(Matrix(A), Matrix(B), Matrix(C))
        ## line below goes with algorithm 1
        # Psi = mysyl(RA,QA,RB,QB,C)
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
        OldPsi = Psi
        ## Algorithm 4 of Bean, O'Reilly, Taylor, 2008, 
        ## Algorithms for the Laplace–Stieltjes Transforms of First Return Times for Stochastic Fluid Flows, 
        ## Methodol Comput Appl Probab (2008) 10:381–408, DOI 10.1007/s11009-008-9077-3
        A = EvalD["+","+"] + Psi * EvalD["-","+"]
        B = EvalD["-","-"] + EvalD["-","+"] * Psi
        C = EvalD["+","-"] - Psi * EvalD["-","+"] * Psi
        ## Algorithm 1 of the above citation
        # A, B dont change, need only comput schur decomp once 
        # C = EvalD["+-"] + Psi * EvalD["-+"] * Psi
    end
    if flag == 1
        exitflag = string(
            "Reached Max Iters ",
            MaxIters,
            " with error ",
            string(maximum(abs.(OldPsi - Psi))),
        )
    end
    v && println("UPDATE: Iterations for Ψ(s=", s,") exited with flag: ", exitflag)
    return Psi
end

function mysyl(RA,QA,RB,QB,C)
    # sylvester but we give it the schur decomps
    # RA, QA = schur(A)
    # RB, QB = schur(B)

    D = -(adjoint(QA) * (C*QB))
    Y, scale = LinearAlgebra.LAPACK.trsyl!('N','N', RA, RB, D)
    LinearAlgebra.rmul!(QA*(Y * adjoint(QB)), inv(scale))
end


"""
Returns the DG approximation to the return probabilities ``ξ`` for the process
``Y(t)``.
NOTE: IMPLEMENTED FOR LAGRANGE BASIS ONLY

    MakeXi(
        B::Generator,
        Ψ::Array;
        mesh=1,
        model=1,
    )

# Arguments
- `B`: an object as returned by `MakeB`
- `Ψ::Array{Float64,2}`: an array as returned by `PsiFun`

# Output
- `ξ::Array{Float64,2}`: a row-vector of first return probabilities
"""
function MakeXi(
    B::Generator,
    Ψ::Array{Float64,2};
    mesh::Mesh = DGMesh(),
    model::Model = Model(),
)
    # BBullet = [B["--"] B["-0"]; B["0-"] B["00"]]
    # invB = inv(Matrix(Bbullet))

    # BBulletPlus = [B["-+"]; B["0+"]]

    # solve the linear system -[ξ 0] [Bmm Bm0; B₀₋ B00]^-1 [B₋₊; B₀₊]Ψ = ξ
    # writing out the system it turns out we only need the index -- and -0
    # blocks of the inverse. Wikipedia tells us that these are

    # the next two lines are a fancy way to do an inverse of a block matrix
    tempMat = inv(Matrix(B[("0","0"),(:,:)]))
    invBmm = inv(B[("-","-"),(:,:)] - B[("-","0"),(:,:)]*tempMat*B[("0","-"),(:,:)])
    invBm0 = -invBmm*B[("-","0"),(:,:)]*tempMat

    A = -(invBmm*B[("-","+"),(:,:)]*Ψ + invBm0*B[("0","+"),(:,:)]*Ψ + LinearAlgebra.I)
    b = zeros(1,size(B[("-","-"),(:,:)],1))

    A[:,1] .= 1.0 # normalisation conditions
    b[1] = 1.0 # normalisation conditions

    ξ = b/A

    return ξ
end

"""
Returns the DG approximation to some quantities regarding the limiting
distribution of a  See Ouput below
NOTE: IMPLEMENTED FOR LAGRANGE BASIS ONLY

    MakeLimitDistMatrices(
        B::Generator,
        D::Dict{Tuple{String,String},Any},
        R::Dict{Tuple{String,Union{Int,Colon}},SparseArrays.SparseMatrixCSC{Float64,Int64}},
        Ψ::Array{<:Real},
        ξ::Array{<:Real},
        mesh::Mesh,
        model::Model;
    )

# Arguments
- `B`: an object as returned by `MakeB`
- `D`: an object as returned by `MakeD`
- `R`: an object as returned by `MakeR`
- `Ψ::Array{Float64,2}`: an array as returned by `PsiFun`
- `ξ::Array{Float64,2}`: an row-vector as returned by `XiFun`

# Output
marginalX, p, K
- `marginalX::Array{Float64,2}`: a row-vector of the marginal limiting
    distribution of the first buffer ``X(t)``.
- `p::Array{Float64,2}`: a row-vector of the distribution of ``X(t)`` at the
    times of first return of ``Y(t)`` to 0.
- `K::Array{Float64,2}`: the array in the operator exponential of the stationary
    distribution

"""
function MakeLimitDistMatrices(
    B::Generator,
    D::Dict{Tuple{String,String},Any},
    R::Dict{Tuple{String,Union{Int,Colon}},SparseArrays.SparseMatrixCSC{Float64,Int64}},
    Ψ::Array{<:Real},
    ξ::Array{<:Real},
    mesh::Mesh;
)
    B00inv = inv(Matrix(B[("0","0"),(:,:)]))
    invBmm = inv(B[("-","-"),(:,:)] - B[("-","0"),(:,:)]*B00inv*B[("0","-"),(:,:)])
    invBm0 = -invBmm*B[("-","0"),(:,:)]*B00inv

    # B00inv = inv(Matrix(B["00"]))
    # invBmm = inv(B["--"] - B["-0"]*B00inv*B["0-"])
    # invBm0 = -invBmm*B["-0"]*B00inv

    αp = ξ * -[invBmm invBm0]

    K = D["+","+"]() + Ψ * D["-","+"]()

    n₊, n₀ = size(B[("+","0"),(:,:)])
    n₋ = size(B[("-","+"),(:,:)],1)

    BBulletPlus = [B[("-","+"),(:,:)]; B[("0","+"),(:,:)]]

    αintegralPibullet = ((αp * BBulletPlus) / -K) * [R["+",:] Ψ*R["-",:]]
    αintegralPi0 = -αintegralPibullet * [B[("+","0"),(:,:)]; B[("-","0"),(:,:)]] * B00inv

    α = sum(αintegralPibullet) + sum(αintegralPi0) + sum(αp)

    p = αp ./ α
    integralPibullet = αintegralPibullet ./ α
    integralPi0 = αintegralPi0 ./ α

    marginalX = zeros(Float64, n₊ + n₋ + n₀)
    idx₊ = [mesh.Fil["p+",:]; repeat(mesh.Fil["+",:]', NBases(mesh), 1)[:]; mesh.Fil["q+",:]]
    marginalX[idx₊] = integralPibullet[1:n₊]
    idx₋ = [mesh.Fil["p-",:]; repeat(mesh.Fil["-",:]', NBases(mesh), 1)[:]; mesh.Fil["q-",:]]
    marginalX[idx₋] = integralPibullet[(n₊+1):end] + p[1:n₋]
    idx₀ = [mesh.Fil["p0",:]; repeat(mesh.Fil["0",:]', NBases(mesh), 1)[:]; mesh.Fil["q0",:]]
    marginalX[idx₀] = integralPi0[:] + p[(n₋+1):end]

    return marginalX, p, K
end
