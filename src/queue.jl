struct Rates
    interior::Matrix{Float64}
    lwr::Vector{Float64}
    upr::Vector{Float64}
end

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