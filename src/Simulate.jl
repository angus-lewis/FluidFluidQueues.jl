"""
Simulates a SFFM defined by `model` until the `StoppingTime` has occured,
given the `InitialCondition` on (φ(0),X(0),Y(0)).

    SimSFFM(
        model::Model,
        StoppingTime::Function,
        InitCondition::NamedTuple{(:φ, :X, :Y)},
    )

# Arguments
- `model`: A Model object
- `StoppingTime`: A function which takes the value of the process at the current
    time and at the time of the last jump of the phase process, as well as the
    `model` object.
    i.e. `StoppingTime(;model,SFFM,SFFM0)` where `SFFM` and `SFFM0` are tuples
    with keys `(:t::Float64, :φ::Int, :X::Float64, :Y::Float64, :n::Int)` which
    are the value of the SFFM at the current time, and time of the previous jump
    of the phase process, repsectively. The `StoppingTime` must return a
    `NamedTuple{(:Ind, :SFFM)}` type where `:Ind` is a `:Bool` value stating
    whether the stopping time has occured or not and `:SFFM` is a tuple in the
    same form as the input `SFFM` but which contains the value of the SFFM at
    the stopping time.
- `InitCondition`: `NamedTuple` with keys `(:φ, :X, :Y)`, `InitCondition.φ` is a
    vector of length `M` of initial states for the phase, `InitCondition.X` is a
    vector of length `M` of initial states for the X-level, `InitCondition.Y` is
    a vector of length `M` of initial states for the Y-level. `M` is the number
    of simulations to be done.

# Output
- a tuple with keys
    - `t::Array{Float64,1}` a vector of length `M` containing the values of
        `t` at the `StoppingTime`.
    - `φ::Array{Float64,1}` a vector of length `M` containing the values of
        `φ` at the `StoppingTime`.
    - `X::Array{Float64,1}` a vector of length `M` containing the values of
        `X` at the `StoppingTime`.
    - `Y::Array{Float64,1}` a vector of length `M` containing the values of
        `Y` at the `StoppingTime`.
    - `n::Array{Float64,1}` a vector of length `M` containing the number of
        transitions of `φ` at the `StoppingTime`
"""
function DiscretisedFluidQueues.simulate(
    ffq::FluidFluidQueue{<:DiscretisedFluidQueue},
    StoppingTime::Function,
    InitCondition::NamedTuple{(:φ, :X, :Y)},
    rng::Random.AbstractRNG=Random.default_rng(),
)
    model = ffq.dq.model
    # find transition probabilities of phase process
    d = LinearAlgebra.diag(model.T)
    P = (model.T - LinearAlgebra.diagm(0 => d)) ./ -d
    # convenient to sample transitions with
    CumP = cumsum(P, dims = 2)

    # transition probabilities upon hitting a boundary
    # lower
    idx₋ = findall(rates(model).<0.0)
    CumP_lwr = zeros(size(P))
    CumP_lwr[idx₋,:] = cumsum(model.P_lwr, dims = 2)
    # upper 
    idx₊ = findall(rates(model).>0.0)
    CumP_upr = zeros(size(P))
    CumP_upr[idx₊,:] = cumsum(model.P_upr, dims = 2)

    # rates for exponential holding times in each phase 
    Λ = LinearAlgebra.diag(model.T)

    # containers to store the sims 
    M = length(InitCondition.φ)
    tSims = Array{Float64,1}(undef, M)
    φSims = Array{Float64,1}(undef, M)
    XSims = Array{Float64,1}(undef, M)
    YSims = Array{Float64,1}(undef, M)
    nSims = Array{Float64,1}(undef, M)

    # simulate M realisations 
    for m = 1:M
        SFM0 = (t=0.0, φ=InitCondition.φ[m], X=InitCondition.X[m], n=0.0)
        Y0 = InitCondition.Y[m]
        # check if initial condition is invalid and return NaNs if it is 
        if !(0.0 <= SFM0.X <= model.b) || !(0.0 <= Y0) || !in(SFM0.φ, 1:n_phases(model))
            (tSims[m], φSims[m], XSims[m], nSims[m]) =
                (t = NaN, φ = NaN, X = NaN, n = NaN)
            YSims[m] = NaN
        else
            cell_idx0 = initialise_cell_idx(ffq,SFM0)
            while true
                # generate random exponential holding time via inverse CDF method
                S = log(rand(rng)) / Λ[SFM0.φ]
                X, φ, t = DiscretisedFluidQueues.Update(model, SFM0, S, CumP, CumP_lwr, CumP_upr, rng)
                cell_idx = find_cell_idx(ffq,SFM0.φ,cell_idx0,X) # find the cell in which X resides
                SFM = (t=t, φ=φ, X=X, n=SFM0.n+1)
                Y = UpdateYt(ffq, SFM, SFM0, Y0, cell_idx, cell_idx0)
                τ = StoppingTime(ffq, SFM, SFM0, Y, Y0)
                if τ.Ind
                    (tSims[m], φSims[m], XSims[m], nSims[m]) = τ.SFM
                    YSims[m] = τ.Y
                    break
                end
                cell_idx0 = cell_idx
                SFM0 = SFM
                Y0 = Y
            end
        end
    end
    return (t = tSims, φ = φSims, X = XSims, n = nSims), YSims
end

function initialise_cell_idx(ffq,SFM0)
    cell_idx0 = findfirst(ffq.dq.mesh.nodes .> SFM0.X)
    if (cell_idx0===nothing)
        if (rates(ffq.dq,SFM0.φ)<0.0) 
            cell_idx0 = length(ffq.dq.mesh.nodes)-1
        else
            cell_idx0 = length(ffq.dq.mesh.nodes)
        end
    elseif (rates(ffq.dq,SFM0.φ)<=0.0)&&(SFM0.X==0.0)
        cell_idx0 = 0
    else 
        cell_idx0 -= 1
    end
    return cell_idx0
end

function find_cell_idx(ffq,φ0,cell_idx0,X)
    # cell_idx = -1
    if (X==0.0)
        cell_idx = 0
    elseif X==ffq.dq.model.b
        cell_idx = length(ffq.dq.mesh.nodes)
    elseif (rates(ffq.dq,φ0)>0.0)
        cell_idx = NaN
        for n in max(1,cell_idx0):length(ffq.dq.mesh.nodes)
            if ffq.dq.mesh.nodes[n]>X
                cell_idx = n-1
                break
            end
        end
    elseif (rates(ffq.dq,φ0)<0.0)
        cell_idx = NaN
        for n in min(length(ffq.dq.mesh.nodes)-1,cell_idx0):-1:1
            if ffq.dq.mesh.nodes[n]<X
                cell_idx = n
                break
            end
        end
    else
        cell_idx = cell_idx0
    end
    return cell_idx
end

"""
Returns ``Y(t+S)`` given ``Y(t)``.

    UpdateYt(
        model::Model,
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
        S::Real,
    )

# Arguments
- `model`: a Model object
- `SFFM0::NamedTuple` containing at least the keys `:X` giving the value of
    ``X(t)`` at the current time, and `:Y` giving the value of ``Y(t)`` at the
    current time, and `:φ` giving the value of `φ(t)`` at the current time.
- `S::Real`: an elapsed amount of time to evaluate ``X`` at, i.e. ``X(t+S)``.
"""
function UpdateYt(
    ffq::FluidFluidQueue,
    SFM::NamedTuple{(:t, :φ, :X, :n)},
    SFM0::NamedTuple{(:t, :φ, :X, :n)},
    Y0::Float64,
    cell_idx::Int, cell_idx0::Int,
)
    # given the position of a SFM, SFM0 at time t, a time step of size s, find the
    # position of Y at time t+s
    speedX = rates(ffq.dq,SFM0.φ)
    if (speedX==0.0)||((speedX<=0.0)&&(SFM0.X==0.0))||((speedX>=0.0)&&(SFM0.X==ffq.dq.model.b))
        # determine if X is at a boundary
        Y = Y0 + (SFM.t-SFM0.t) * getrates(ffq,SFM0.φ,cell_idx0)
    else
        # determine if X has hit a boundary (did not start at the same boundary)
        if (SFM.X==0.0)||(SFM.X==ffq.dq.model.b) 
            if speedX>0.0 # must be upper boundary in this case
                cells_traversed_fully = cell_idx0+1:cell_idx-1
                time_spent_in_cells_traversed = [
                    (ffq.dq.mesh.nodes[cell_idx0+1]-SFM0.X)./speedX; # time to leave cell_idx0
                    [Δ(ffq.dq.mesh,k) for k in cells_traversed_fully]./speedX;
                    #(SFM.X-ffq.dq.mesh.nodes[cell_idx])./speedX; # time since last change of cell
                ]
                Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx0:cell_idx-1),time_spent_in_cells_traversed)
            elseif speedX<0.0 # must be lower boundary in this case
                cells_traversed_fully = cell_idx+1:cell_idx0-1
                time_spent_in_cells_traversed = -[
                    #(ffq.dq.mesh.nodes[cell_idx+1]-SFM.X)./speedX; # time since last change of cell
                    [Δ(ffq.dq.mesh,k) for k in cells_traversed_fully]./speedX;
                    (SFM0.X-ffq.dq.mesh.nodes[cell_idx0])./speedX; # time to leave first cell
                ]
                Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx+1:cell_idx0),time_spent_in_cells_traversed)
            end
        else
            if speedX>0.0
                if cell_idx0==cell_idx
                    time_spent_in_cells_traversed = (SFM.X-SFM0.X)./speedX
                    Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx0),time_spent_in_cells_traversed)
                else
                    cells_traversed_fully = cell_idx0+1:cell_idx-1
                    time_spent_in_cells_traversed = [
                        (ffq.dq.mesh.nodes[cell_idx0+1]-SFM0.X)./speedX; # time to leave cell_idx0
                        [Δ(ffq.dq.mesh,k) for k in cells_traversed_fully]./speedX;
                        (SFM.X-ffq.dq.mesh.nodes[cell_idx])./speedX; # time since last change of cell
                    ]
                    Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx0:cell_idx),time_spent_in_cells_traversed)
                end
            elseif speedX<0.0
                if cell_idx0==cell_idx
                    time_spent_in_cells_traversed = -(SFM0.X-SFM.X)./speedX
                    Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx0),time_spent_in_cells_traversed)
                else
                    cells_traversed_fully = cell_idx+1:cell_idx0-1
                    time_spent_in_cells_traversed = -[
                        (ffq.dq.mesh.nodes[cell_idx+1]-SFM.X)./speedX; # time since last change of cell
                        [Δ(ffq.dq.mesh,k) for k in cells_traversed_fully]./speedX;
                        (SFM0.X-ffq.dq.mesh.nodes[cell_idx0])./speedX; # time to leave first cell
                    ]
                    Y = Y0 + dot(getrates(ffq,SFM0.φ,cell_idx:cell_idx0),time_spent_in_cells_traversed)
                end
            end
        end
    end
    return Y
end


"""
Constructs the `StoppingTime` which is the first exit of the process ``Y(t)``
from the interval ``[u,v]``. ASSUMES ``Y(t)`` is monotonic between jumps.

    first_exit_y( u::Real, v::Real)

# Arguments
- `u`: a lower boundary
- `v`: an upper boundary

# Output
- `FirstExitYFun`: a function with one method
    - `FirstExitYFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )`: a stopping time for a SFFM
"""
function first_exit_y( u::Real, v::Real)
    # SFFM Method
    function first_exit_yFun(
        ffq::FluidFluidQueue,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
        Y::Float64,
        Y0::Float64,
    )
        Ind = ((Y < u) || (Y > v))
        if Ind
            idx = [Y < u; Y > v]
            boundaryHit = only([u;v][idx])
            cell_idx0 = initialise_cell_idx(ffq,SFM0)
            SFM_path(t) = (t=SFM0.t+t,φ=SFM0.φ,X=DiscretisedFluidQueues.UpdateX(ffq.dq.model,SFM0,t),n=SFM0.n)
            function YFun(t)
                SFMt = SFM_path(t)
                cell_idx = find_cell_idx(ffq,SFM0.φ,cell_idx0,SFMt.X)
                return UpdateYt(ffq, SFMt, SFM0, Y0, cell_idx, cell_idx0) - boundaryHit
            end
            S = SFM.t - SFM0.t
            tstar = fzero(YFun, 1e-64, S)
            X = DiscretisedFluidQueues.UpdateX(ffq.dq.model, SFM0, tstar)
            t = SFM0.t + tstar
            Y = YFun(tstar)+boundaryHit
            SFM = (t, SFM0.φ, X, SFM0.n)
        end
        return (Ind = Ind, SFM = SFM, Y=Y)
    end
    return first_exit_yFun
end

function n_jumps_y(N::Int)
    # SFFM Method
    function n_jumps_yFun(
        ffq::FluidFluidQueue,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
        Y::Float64,
        Y0::Float64,
    )
        Ind = SFM.n>=N
        return (Ind = Ind, SFM = SFM, Y=Y)
    end
    return n_jumps_yFun
end

function fixed_time_y(T::Real)
    # SFFM Method
    function fixed_time_yFun(
        ffq::FluidFluidQueue,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
        Y::Float64,
        Y0::Float64,
    )
        Ind = SFM.t>=T
        if Ind 
            S = T-SFM0.t
            X = DiscretisedFluidQueues.UpdateX(ffq.dq.model,SFM0,S)
            cell_idx = find_cell_idx(ffq,SFM0.φ,cell_idx0,X) # find the cell in which X resides
            SFM = (t=T, φ=SFM0.φ, X=X, n=SFM0.n)
            Y = UpdateYt(ffq, SFM, SFM0, Y0, cell_idx, cell_idx0)
        end 
        return (Ind = Ind, SFM = SFM, Y=Y)
    end
    return fixed_time_yFun
end


"""
Finds zero of `f` using the bisection method on the interval `[a,b]` with
error `err`.

    fzero( f::Function, a::Real, b::Real; err::Float64 = 1e-8)

"""
function fzero( f::Function, a::Real, b::Real; err::Float64 = 1e-6)
    # finds zeros of f using the bisection method
    c = a + (b - a) / 2.0
    while a < c < b
        fc = f(c)
        if (abs(fc) < err)&&((b-a)<err)
            break
        end
        p = f(a) * fc 
        if p < 0.0
            a, b = a, c
        elseif p > 0.0
            a, b = c, b
        elseif fc == 0.0
            break
        else 
            a = a+sqrt(eps())
        end
        c = a + (b - a) / 2.0
    end
    return c
end

function brute_force_zero(f::Function, a::Real, b::Real; err::Float64 = 1e-6)
    x = range(a,b;step=err)
    fx = abs.(f.(x))
    val, ind = findmin(fx)
    (val>err)&&(@warn "no values with |f(x)| sufficiently small, returning the closest: |f(x)| = f("*string(x[ind])*") = "*string(fx[ind]))
    return x[ind]
end
