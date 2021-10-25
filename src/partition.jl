###############
###############
## partition ##
###############
###############

abstract type AbstractSign end 
struct Plus <: AbstractSign end
struct Minus <: AbstractSign end
struct Zero <: AbstractSign end
struct NotPlus <: AbstractSign end
struct NotMinus <: AbstractSign end
struct NotZero <: AbstractSign end

function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{Plus},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.lwr[bndry_idx]>0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{Minus},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.lwr[bndry_idx]<0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{Zero},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.lwr[bndry_idx]==0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{NotPlus},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.lwr[bndry_idx]>0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{NotMinus},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.lwr[bndry_idx]<0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{NotZero},i::Int)
    if DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.lwr[bndry_idx]==0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end

function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{Plus},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.upr[bndry_idx]>0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{Minus},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.upr[bndry_idx]<0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{Zero},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if ffq.rates.upr[bndry_idx]==0.0
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{NotPlus},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.upr[bndry_idx]>0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{NotMinus},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.upr[bndry_idx]<0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end
function _sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{NotZero},i::Int)
    if DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.generator.dq.model.S[1:i]))
        if !(ffq.rates.upr[bndry_idx]==0.0)
            return true
        else 
            return false
        end
    else
        return false
    end
end

_sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Type{<:AbstractSign}) = 
    [_sign_map_lwr_bndry(ffq,k,i) for i in 1:n_phases(ffq.generator.dq)]
_sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{<:AbstractSign}) = 
    [_sign_map_upr_bndry(ffq,k,i) for i in 1:n_phases(ffq.generator.dq)]
_sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Colon) = 
    [_sign_map_lwr_bndry(ffq,NotPlus,i)||_sign_map_lwr_bndry(ffq,Plus,i) for i in 1:n_phases(ffq.generator.dq)]
_sign_map_upr_bndry(ffq::FluidFluidQueue,k::Colon) = 
    [_sign_map_upr_bndry(ffq,NotPlus,i)||_sign_map_upr_bndry(ffq,Plus,i) for i in 1:n_phases(ffq.generator.dq)]


_sign_map_interior(ffq::FluidFluidQueue,k::Type{Plus}) = map(x->x>0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{Minus}) = map(x->x<0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{Zero}) = map(x->x==0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotPlus}) = map(x->!(x>0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotMinus}) = map(x->!(x<0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotZero}) = map(x->!(x==0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Colon) = trues(size(ffq.rates.interior))

_sign_map(ffq::FluidFluidQueue,k) = 
    [_sign_map_lwr_bndry(ffq,k) _sign_map_interior(ffq,k) _sign_map_upr_bndry(ffq,k)]


Base.size(R::FluidFluidOperator{1}) = size(R.r)
Base.getindex(R::FluidFluidOperator{1},i) = R.r[i]
# function Base.getindex(R::FluidFluidOperator{1},(ks,is,ps)::Tuple{Any,Any,Any})
#     return Base.getindex(R,(ks,is,ps)...)
# end
function Base.getindex(R::FluidFluidOperator{1},(ks,is,ps)::Tuple{Any,Any,Any})

    (typeof(is)==Colon)&&(is=1:n_phases(R.ffq.generator.dq))
    (typeof(ks)==Colon)&&(ks=0:n_intervals(R.ffq.generator.dq)+1)
    (typeof(ps)==Colon)&&(ps=1:n_bases_per_cell(R.ffq.generator.dq))

    rows = Int[]
    for k in ks 
        for i in is
            for p in ps 
                if (k==0)&&(p==1)&&DiscretisedFluidQueues._has_left_boundary(R.ffq.generator.dq.model.S,i)
                    push!(rows,DiscretisedFluidQueues._map_to_index((0,i,1),R.ffq.generator.dq))
                elseif (k==n_intervals(R.ffq.generator.dq)+1)&&(p==1)&&DiscretisedFluidQueues._has_right_boundary(R.ffq.generator.dq.model.S,i)
                    push!(rows,DiscretisedFluidQueues._map_to_index((n_intervals(R.ffq.generator.dq)+1,i,1),R.ffq.generator.dq))
                elseif (k!=0)&&(k!=n_intervals(R.ffq.generator.dq)+1)
                    push!(rows,DiscretisedFluidQueues._map_to_index((k,i,p),R.ffq.generator.dq))
                end
            end
        end
    end
    return R.r[rows]
end

function Base.getindex(R::FluidFluidOperator{1},(k,i,p)::Tuple{DataType,Any,Any})
    k_map = findall(_sign_map(R.ffq,k))
    (i==:)&&(i=1:n_phases(R.ffq.generator.dq))
    (p==:)&&(p=1:n_bases_per_cell(R.ffq.generator.dq))
    k_map = k_map[map(x->x[1]∈i,k_map)]
    return vcat([R[(k_map[row][2]-1,k_map[row][1],p)] for row in 1:length(k_map)]...)
end
Base.getindex(R::FluidFluidOperator{1},(k,i)::Tuple{Any,Any}) = R[(k,i,:)]
Base.getindex(R::FluidFluidOperator{1},k::Type{<:AbstractSign}) = R[(k,:,:)]

Base.size(ffq::FluidFluidOperator{2}) = size(ffq.generator)
Base.getindex(ffq::FluidFluidOperator{2},i,j) = ffq.generator[i,j]

function Base.getindex(ffq::FluidFluidOperator{2},(k,i,p)::Tuple{Any,Any,Any},(l,j,q)::Tuple{Any,Any,Any})
    k_map = findall(_sign_map(ffq,k))
    l_map = findall(_sign_map(ffq,l))
    (i==:)&&(i=1:n_phases(ffq.generator.dq))
    (j==:)&&(j=1:n_phases(ffq.generator.dq))
    (p==:)&&(p=1:n_bases_per_cell(ffq.generator.dq))
    (q==:)&&(q=1:n_bases_per_cell(ffq.generator.dq))
    k_map = k_map[map(x->x[1]∈i,k_map)]
    l_map = l_map[map(x->x[1]∈j,l_map)]
    return vcat([hcat([ffq.generator[(k_map[row][2]-1,k_map[row][1],p),(l_map[col][2]-1,l_map[col][1],q)] for col in 1:length(l_map)]...) for row in 1:length(k_map)]...)
end

Base.getindex(ffq::FluidFluidOperator{2},(k,i)::Tuple{Any,Any},(l,j)::Tuple{Any,Any}) = ffq[(k,i,:),(l,j,:)]
Base.getindex(ffq::FluidFluidOperator{2},k::Type{<:AbstractSign},l) = ffq[(k,:,:),(l,:,:)]
