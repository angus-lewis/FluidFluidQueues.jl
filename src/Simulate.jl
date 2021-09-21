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
function SimSFFM(
    model::Model,
    StoppingTime::Function,
    InitCondition::NamedTuple{(:φ, :X, :Y)},
)
    d = LinearAlgebra.diag(model.T)
    P = (model.T - LinearAlgebra.diagm(0 => d)) ./ -d
    CumP = cumsum(P, dims = 2)
    Λ = LinearAlgebra.diag(model.T)

    M = length(InitCondition.φ)
    tSims = Array{Float64,1}(undef, M)
    φSims = Array{Float64,1}(undef, M)
    XSims = Array{Float64,1}(undef, M)
    YSims = Array{Float64,1}(undef, M)
    nSims = Array{Float64,1}(undef, M)

    for m = 1:M
        SFFM0 = (
            t = 0.0,
            φ = InitCondition.φ[m],
            X = InitCondition.X[m],
            Y = InitCondition.Y[m],
            n = 0.0,
        )
        if !(model.Bounds[1, 1] <= SFFM0.X <= model.Bounds[1, 2]) ||
           !(model.Bounds[2, 1] <= SFFM0.Y <= model.Bounds[2, 2]) ||
           !in(SFFM0.φ, 1:n_phases(model))
            (tSims[m], φSims[m], XSims[m], YSims[m], nSims[m]) =
                (t = NaN, φ = NaN, X = NaN, Y = NaN, n = NaN)
        else
            while 1 == 1
                S = log(rand()) / Λ[SFFM0.φ]
                t = SFFM0.t + S
                X = UpdateXt(model, SFFM0, S)
                Y = UpdateYt(model, SFFM0, S)
                φ = findfirst(rand() .< CumP[SFFM0.φ, :])
                n = SFFM0.n + 1.0
                SFFM = (t = t, φ = φ, X = X, Y = Y, n = n)
                τ = StoppingTime(model, SFFM, SFFM0)
                if τ.Ind
                    (tSims[m], φSims[m], XSims[m], YSims[m], nSims[m]) = τ.SFFM
                    break
                end
                SFFM0 = SFFM
            end
        end
    end
    return (t = tSims, φ = φSims, X = XSims, Y = YSims, n = nSims)
end

"""
Returns ``Y(t+S)`` given ``Y(t)``.

    UpdateYt(
        model::Model,
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
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
    model::Model,
    SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    S::Real,
)
    # given the last position of a SFFM, SFFM0, a time step of size s, find the
    # position of Y at time t
    if rates(model,SFFM0.φ) == 0
        Y = SFFM0.Y + S * model.r.r(SFFM0.X)[SFFM0.φ]
    else
        X = UpdateXt(model, SFFM0, S)
        ind = (X.==model.Bounds[1, :])[:]
        if any(ind)
            # time at which Xt hits u or v
            t0 = (model.Bounds[1, ind][1] - SFFM0.X) / rates(model,SFFM0.φ)
            Y =
                SFFM0.Y +
                (model.r.R(X)[SFFM0.φ] - model.r.R(SFFM0.X)[SFFM0.φ]) /
                rates(model,SFFM0.φ) +
                model.r.r(model.Bounds[1, ind][1])[SFFM0.φ] * (S - t0)
        else
            Y =
                SFFM0.Y +
                (model.r.R(X)[SFFM0.φ] - model.r.R(SFFM0.X)[SFFM0.φ]) /
                rates(model,SFFM0.φ)
        end
    end
    return Y
end


"""
Constructs the `StoppingTime` ``1(t>T)``

    FixedTime( T::Real)

# Arguments
- `T`: a time at which to stop the process

# Output
- `FixedTimeFun`: a function with two methods
    - `FixedTimeFun(
        model::Model,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
    )`: a stopping time for a SFM.
    - `FixedTimeFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )`: a stopping time for a SFFM
"""
function FixedTime( T::Float64)
    # Defines a simple stopping time, 1(t>T).
    # SFFM METHOD
    function FixedTimeFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )
        Ind = SFFM.t > T
        if Ind
            s = T - SFFM0.t
            X = UpdateXt(model, SFFM0, s)
            Y = UpdateYt(model, SFFM0, s)
            SFFM = (T, SFFM0.φ, X, Y, SFFM0.n)
        end
        return (Ind = Ind, SFFM = SFFM)
    end
    return FixedTimeFun
end


"""
Constructs the `StoppingTime` ``1(N(t)>n)`` where ``N(t)`` is the number of
jumps of ``φ`` by time ``t``.

    NJumps( N::Int)

# Arguments
- `N`: a desired number of jumps

# Output
- `NJumpsFun`: a function with two methods
    - `NJumpsFun(
        model::Model,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
    )`: a stopping time for a SFM.
    - `NJumpsFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )`: a stopping time for a SFFM
"""
function NJumps( N::Int)
    # Defines a simple stopping time, 1(n>N), where n is the number of jumps of φ.
    # SFFM method
    function NJumpsFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )
        Ind = n >= N
        return (Ind = Ind, SFFM = SFFM)
    end
    return NJumpsFun
end

"""
Constructs the `StoppingTime` which is the first exit of the process ``X(t)``
from the interval ``[u,v]``.

    FirstExitX( u::Real, v::Real)

# Arguments
- `u`: a lower boundary
- `v`: an upper boundary

# Output
- `FirstExitXFun`: a function with two methods
    - `FirstExitXFun(
        model::Model,
        SFM::NamedTuple{(:t, :φ, :X, :n)},
        SFM0::NamedTuple{(:t, :φ, :X, :n)},
    )`: a stopping time for a SFM.
    - `FirstExitXFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )`: a stopping time for a SFFM
"""
function FirstExitX( u::Real, v::Real)
    # SFFM Method
    function FirstExitXFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )
        Ind = ((X > v) || (X < u))
        if Ind
            if (X > v)
                X = v
            else
                X = u
            end
            s = (X - SFFM0.X) / rates(model,SFFM0.φ) # can't exit with c = 0.
            t = SFFM0.t + s
            Y = UpdateYt(model, SFFM0, s)
            SFFM = (t, SFFM0.φ, X, Y, SFFM0.n)
        end
        return (Ind = Ind, SFFM = SFFM)
    end
    return FirstExitXFun
end

"""
Constructs the `StoppingTime` which is the first exit of the process ``Y(t)``
from the interval ``[u,v]``. ASSUMES ``Y(t)`` is monotonic between jumps.

    FirstExitY( u::Real, v::Real)

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
function FirstExitY( u::Real, v::Real)
    # SFFM Method
    function FirstExitYFun(
        model::Model,
        SFFM::NamedTuple{(:t, :φ, :X, :Y, :n)},
        SFFM0::NamedTuple{(:t, :φ, :X, :Y, :n)},
    )
        Ind = ((Y < u) || (Y > v))
        if Ind
            idx = [Y < u; Y > v]
            boundaryHit = [u;v][idx][1]
            YFun(t) = UpdateYt(model, SFFM0, t) - boundaryHit
            S = t - SFFM0.t
            tstar = fzero(YFun, 0, S)
            X = UpdateXt(model, SFFM0, tstar)
            t = SFFM0.t + tstar
            Y = boundaryHit
            SFFM = (t, SFFM0.φ, X, Y, SFFM0.n)
        end
        return (Ind = Ind, SFFM = SFFM)
    end
    return FirstExitYFun
end


"""
Finds zero of `f` using the bisection method on the interval `[a,b]` with
error `err`.

    fzero( f::Function, a::Real, b::Real; err::Float64 = 1e-8)

"""
function fzero( f::Function, a::Real, b::Real; err::Float64 = 1e-6)
    # finds zeros of f using the bisection method
    c = a + (b - a) / 2
    while a < c < b
        fc = f(c)
        if abs(fc) < err
            break
        end
        p = f(a) * fc 
        if p < 0
            a, b = a, c
        elseif p > 0
            a, b = c, b
        elseif fc == 0
            break
        else 
            a = a+sqrt(eps())
        end
        c = a + (b - a) / 2
    end
    return c
end
