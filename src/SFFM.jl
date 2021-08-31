struct FluidFluidQueue <: Model
    T::Array{<:Real,2}
    S::PhaseSet
    r::NamedTuple{(:r, :R, :a)}
    Bounds::Array{<:Real}
end 
get_rates(m::FluidFluidQueue) = get_rates(m.S)
get_rates(m::FluidFluidQueue,i::Int) = get_rates(m.S,i)
n_phases(m::FluidFluidQueue) = n_phases(m.S)
phases(m::FluidFluidQueue) = 1:n_phases(m.S)

function FluidFluidQueue(
    T::Array{<:Real},
    C::Array{<:Real,1},
    r::NamedTuple{(:r, :R)};
    Bounds::Array{<:Real,2} = [-Inf Inf; -Inf Inf],
    v::Bool = false,
)
    a(x) = abs.(r.r(x))
    r = (r = r.r, R = r.R, a = a)

    v && println("UPDATE: FluidFluidQueue object created with fields ", fieldnames(Model))
    return FluidFluidQueue(T,PhaseSet(C),r,Bounds)
end
function FluidFluidQueue(
    T::Array{<:Real},
    S::PhaseSet,
    r::NamedTuple{(:r, :R)};
    Bounds::Array{<:Real,2} = [-Inf Inf; -Inf Inf],
    v::Bool = false,
)
    a(x) = abs.(r.r(x))
    r = (r = r.r, R = r.R, a = a)

    v && println("UPDATE: FluidFluidQueue object created with fields ", fieldnames(Model))
    return FluidFluidQueue(T,S,r,Bounds)
end
FluidFluidQueue() = FluidFluidQueue([0],PhaseSet([0]),(r=0, R=0, a=0),[0])

function _duplicate_zero_states(T::Array{<:Real,2},C::Array{<:Real,1}, r::NamedTuple{(:r, :R)})
    T_aug, C_aug = _duplicate_zero_states(T,C)
    
    plus_idx = falses(length(C_aug)) # zero states associated with +
    neg_idx = falses(length(C_aug)) # zero states associated with -
    for i in 1:length(C)
        if C[i] == 0
            plus_idx[i+c_zero] = true
            neg_idx[i+c_zero+1] = true
            c_zero += 1
        end
    end

    # assign second fluid rates
    function r_aug_inner(x)
        out = zeros(length(C_aug))
        out[(C_aug.!=0).|(plus_idx)] = r.r(x)
        out[neg_idx] = r.r(x)[C.==0]
        return out
    end
    function R_aug(x)
        out = zeros(length(C_aug))
        out[(C_aug.!=0).|(plus_idx)] = r.R(x)
        out[neg_idx] = r.R(x)[C.==0]
        return out
    end
    r_aug = (r = r_aug_inner, R = R_aug)

    return T_aug, C_aug, r_aug
end

function augment_model(model::FluidFluidQueue)
    if (any(model.C.==0))
        T_aug, C_aug, r_aug = _duplicate_zero_states(model.T,model.C,model.r)
        return FluidFluidQueue(T_aug,C_aug,r_aug,model.Bounds)
    else # no zero states, no augmentation needed
       return model 
    end
end