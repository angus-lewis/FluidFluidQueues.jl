struct Rates 
    interior::Matrix{Float64}
    lwr::Vector{Float64}
    upr::Vector{Float64}
end

"""
Indexing for Rates by cell, k, then phase, i.
"""
function Base.getindex(s::Rates, k::Int, i::Int)
    if k==0
        r = s.lwr[i]
    elseif k==size(s.interior,2)+1
        r = s.upr[i]
    else
        r = s.interior[i,k]
    end
    return r
end
Base.getindex(s::Rates, i::AbstractVector{Int}, k::AbstractVector{Int}) = 
    [s[ii,kk] for ii in i, kk in k]
Base.getindex(s::Rates, i::AbstractVector{Int}, k::Int) = 
    [s[ii,k] for ii in i]
Base.getindex(s::Rates, i::Int, k::AbstractVector{Int}) = 
    [s[i,kk] for kk in k]

struct FluidFluidQueue{T<:DiscretisedFluidQueue}# where N<:AbstractVector{Float64}}
    dq::T
    rates::Rates
    function FluidFluidQueue{T}(dq::T,rates::Rates) where T<:DiscretisedFluidQueue
        if !(n_intervals(dq)==size(rates.interior,2))
            throw(DomainError("size(rates.interior,2) needs to match the number of intervals in the DiscretisedFluidQueue"))
        elseif !(n_phases(dq)==size(rates.interior,1))
            throw(DomainError("size(rates.interior,1) needs to match the number of phases in the DiscretisedFluidQueue"))
        end
        if !(sum(DiscretisedFluidQueues._has_left_boundary(dq.model.S))==length(rates.lwr))
            throw(DomainError("rates.lwr must have length equal to the number of phases with 
                the a left point mass, sum(DiscretisedFluidQueues._has_left_boundary(B.dq.model.S))"))
        end
        if !(sum(DiscretisedFluidQueues._has_right_boundary(dq.model.S))==length(rates.upr))
            throw(DomainError("rates.upr must be have length equal to the number of phases with 
                the a left point mass, sum(DiscretisedFluidQueues._has_left_boundary(B.dq.model.S))"))
        end
        return new{T}(dq,rates)
    end
end 
FluidFluidQueue(dq::DiscretisedFluidQueue,rates::Rates)=FluidFluidQueue{typeof(dq)}(dq,rates)

function FluidFluidQueue(dq::DiscretisedFluidQueue,interior::Matrix{Float64},lwr::Vector{Float64},upr::Vector{Float64})
    return FluidFluidQueue(dq,Rates(interior,lwr,upr))
end

function getrates(ffq::FluidFluidQueue,i::Int,k::Int)
    if k==0
        if rates(ffq.dq,i)<=0.0
            idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
            r = ffq.rates.lwr[idx]
        else
            r = 0.0
        end
    elseif k==size(ffq.rates.interior,2)+1
        if rates(ffq.dq,i)<=0.0
            idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
            r = ffq.rates.upr[idx]
        else
            r = 0.0
        end
    else
        r = ffq.rates.interior[i,k]
    end
    return r
end
getrates(ffq::FluidFluidQueue,i::AbstractVector{Int},k::AbstractVector{Int}) = 
    [getrates(ffq::FluidFluidQueue,ii,kk) for ii in i, kk in  k]
getrates(ffq::FluidFluidQueue,i::Int,k::AbstractVector{Int}) = 
    [getrates(ffq::FluidFluidQueue,i,kk) for kk in  k]
getrates(ffq::FluidFluidQueue,i::AbstractVector{Int},k::Int) = 
    [getrates(ffq::FluidFluidQueue,ii,k) for ii in i]

function DiscretisedFluidQueues.augment_model(ffq::FluidFluidQueue)
    aug_model = augment_model(ffq.dq.model)
    # aug_rates = zeros(n_phases(augment_model),n_intervals(ffq.dq))
    interior_rates = zeros(0,n_intervals(ffq.dq))
    for i in 1:n_phases(ffq.dq.model)
        if rates(ffq.dq)[i]==0.0
            interior_rates = [interior_rates; repeat(ffq.rates.interior[[i],:],2)]
        else
            interior_rates = [interior_rates; ffq.rates.interior[[i],:]]
        end
    end
    aug_dq = DiscretisedFluidQueue(aug_model,ffq.dq.mesh)
    return FluidFluidQueue(aug_dq,interior_rates,ffq.rates.lwr,ffq.rates.upr)
end