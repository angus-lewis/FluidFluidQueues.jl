abstract type FluidFluidOperator{N} <: AbstractArray{Float64,N} end

struct Rates
    interior::Matrix{Float64}
    lwr::Vector{Float64}
    upr::Vector{Float64}
end

struct FluidFluidQueue <: FluidFluidOperator{2}
    generator::Generator
    rates::Rates
    function FluidFluidQueue(B::Generator,rates::Rates)
        if !(n_intervals(B.dq)==size(rates.interior,2))
            throw(DomainError("size(rates.interior,2) needs to match the number of intervals in the DiscretisedFluidQueue"))
        elseif !(n_phases(B.dq)==size(rates.interior,1))
            throw(DomainError("size(rates.interior,1) needs to match the number of phases in the DiscretisedFluidQueue"))
        end
        if !(sum(DiscretisedFluidQueues._has_left_boundary(B.dq.model.S))==length(rates.lwr))
            throw(DomainError("rates.lwr must have length equal to the number of phases with 
                the a left point mass, sum(DiscretisedFluidQueues._has_left_boundary(B.dq.model.S))"))
        end
        if !(sum(DiscretisedFluidQueues._has_right_boundary(B.dq.model.S))==length(rates.upr))
            throw(DomainError("rates.upr must be have length equal to the number of phases with 
                the a left point mass, sum(DiscretisedFluidQueues._has_left_boundary(B.dq.model.S))"))
        end
        return new(B,rates)
    end
end 

function FluidFluidQueue(B::Generator,interior::Matrix{Float64},lwr::Vector{Float64},upr::Vector{Float64})
    return FluidFluidQueue(B,Rates(interior,lwr,upr))
end

struct RatesOperator <: FluidFluidOperator{1} 
    r::Vector{Float64}
    ffq::FluidFluidQueue
end

RatesOperator(ffq::FluidFluidQueue) = RatesOperator(
    abs.([ffq.rates.lwr[:]; 
    repeat(ffq.rates.interior[:],inner=n_bases_per_cell(ffq.generator.dq));
    ffq.rates.upr[:]]),
    ffq)
