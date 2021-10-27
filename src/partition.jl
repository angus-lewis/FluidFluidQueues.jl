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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_left_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    if DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S,i)
        bndry_idx = sum(DiscretisedFluidQueues._has_right_boundary(ffq.dq.model.S[1:i]))
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
    [_sign_map_lwr_bndry(ffq,k,i) for i in 1:n_phases(ffq.dq)]
_sign_map_upr_bndry(ffq::FluidFluidQueue,k::Type{<:AbstractSign}) = 
    [_sign_map_upr_bndry(ffq,k,i) for i in 1:n_phases(ffq.dq)]
_sign_map_lwr_bndry(ffq::FluidFluidQueue,k::Colon) = 
    [_sign_map_lwr_bndry(ffq,NotPlus,i)||_sign_map_lwr_bndry(ffq,Plus,i) for i in 1:n_phases(ffq.dq)]
_sign_map_upr_bndry(ffq::FluidFluidQueue,k::Colon) = 
    [_sign_map_upr_bndry(ffq,NotPlus,i)||_sign_map_upr_bndry(ffq,Plus,i) for i in 1:n_phases(ffq.dq)]


_sign_map_interior(ffq::FluidFluidQueue,k::Type{Plus}) = map(x->x>0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{Minus}) = map(x->x<0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{Zero}) = map(x->x==0.0,ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotPlus}) = map(x->!(x>0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotMinus}) = map(x->!(x<0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Type{NotZero}) = map(x->!(x==0.0),ffq.rates.interior)
_sign_map_interior(ffq::FluidFluidQueue,k::Colon) = trues(size(ffq.rates.interior))

_sign_map(ffq::FluidFluidQueue,k) = 
    [_sign_map_lwr_bndry(ffq,k) _sign_map_interior(ffq,k) _sign_map_upr_bndry(ffq,k)]


Base.size(R::FluidFluidOperator{1}) = size(R.array)
Base.getindex(R::FluidFluidOperator{1},i::AbstractVector) = R.array[i]
Base.getindex(R::FluidFluidOperator{1},i::Int) = R.array[i]
# function Base.getindex(R::FluidFluidOperator{1},(ks,is,ps)::Tuple{Any,Any,Any})
#     return Base.getindex(R,(ks,is,ps)...)
# end
function Base.getindex(R::FluidFluidOperator{1},(ks,is,ps)::Tuple{Any,Any,Any})
    rows = DiscretisedFluidQueues._build_index_from_tuple(R.ffq.dq,(ks,is,ps))
    # (typeof(is)==Colon)&&(is=1:n_phases(R.ffq.dq))
    # (typeof(ks)==Colon)&&(ks=0:n_intervals(R.ffq.dq)+1)
    # (typeof(ps)==Colon)&&(ps=1:n_bases_per_cell(R.ffq.dq))

    # rows = Int[]
    # for k in ks 
    #     for i in is
    #         for p in ps 
    #             if (k==0)&&(p==1)&&DiscretisedFluidQueues._has_left_boundary(R.ffq.dq.model.S,i)
    #                 push!(rows,DiscretisedFluidQueues._map_to_index((0,i,1),R.ffq.dq))
    #             elseif (k==n_intervals(R.ffq.dq)+1)&&(p==1)&&DiscretisedFluidQueues._has_right_boundary(R.ffq.dq.model.S,i)
    #                 push!(rows,DiscretisedFluidQueues._map_to_index((n_intervals(R.ffq.dq)+1,i,1),R.ffq.dq))
    #             elseif (k!=0)&&(k!=n_intervals(R.ffq.dq)+1)
    #                 push!(rows,DiscretisedFluidQueues._map_to_index((k,i,p),R.ffq.dq))
    #             end
    #         end
    #     end
    # end
    return R.array[rows]
end
function _build_index_map_from_tuple_sign(ffq,(k,i,p))
    k_map = findall(_sign_map(ffq,k))
    (i==:)&&(i=1:n_phases(ffq.dq))
    (p==:)&&(p=1:n_bases_per_cell(ffq.dq))
    k_map = k_map[map(x->x[1]∈i,k_map)]
    return k_map
end
function _build_index_from_tuple_sign(ffq,(k,i,p))
    k_map = _build_index_map_from_tuple_sign(ffq,(k,i,p))
    idx = vcat([DiscretisedFluidQueues._build_index_from_tuple(ffq.dq,(k_map[row][2]-1,k_map[row][1],p)) for row in 1:length(k_map)]...)
    return idx 
end
function Base.getindex(R::FluidFluidOperator{1},(k,i,p)::Tuple{DataType,Any,Any})
    idx = _build_index_from_tuple_sign(R.ffq,(k,i,p))
    return R[idx]
end
Base.getindex(R::FluidFluidOperator{1},(k,i)::Tuple{Any,Any}) = R[(k,i,:)]
Base.getindex(R::FluidFluidOperator{1},k::Type{<:AbstractSign}) = R[(k,:,:)]

Base.size(B::FluidFluidOperator{2}) = size(B.array)
Base.getindex(B::FluidFluidOperator{2},i,j) = B.array[i,j]

Base.getindex(B::FluidFluidOperator{2},(k,i,p)::Tuple{Any,Any,Any},(l,j,q)::Tuple{Any,Any,Any}) = 
    DiscretisedFluidQueues._getindex_from_tuple(B,B.ffq.dq,(k,i,p),(l,j,q))
function _getindex_from_tuple_sign(B::FluidFluidOperator{2},(k,i,p)::Tuple,(l,j,q)::Tuple)
    # k_map = findall(_sign_map(B.ffq,k))
    # l_map = findall(_sign_map(B.ffq,l))
    # (i==:)&&(i=1:n_phases(B.ffq.dq))
    # (j==:)&&(j=1:n_phases(B.ffq.dq))
    # (p==:)&&(p=1:n_bases_per_cell(B.ffq.dq))
    # (q==:)&&(q=1:n_bases_per_cell(B.ffq.dq))
    # k_map = k_map[map(x->x[1]∈i,k_map)]
    # l_map = l_map[map(x->x[1]∈j,l_map)]
    # return vcat([hcat([B[(k_map[row][2]-1,k_map[row][1],p),(l_map[col][2]-1,l_map[col][1],q)] for col in 1:length(l_map)]...) for row in 1:length(k_map)]...)
    rows = _build_index_from_tuple_sign(B.ffq,(k,i,p))
    cols = _build_index_from_tuple_sign(B.ffq,(l,j,q))
    return B[rows,cols]
end
Base.getindex(B::FluidFluidOperator{2},(k,i,p)::Tuple{DataType,Any,Any},(l,j,q)::Tuple{Any,Any,Any}) = 
    _getindex_from_tuple_sign(B,(k,i,p),(l,j,q))
Base.getindex(B::FluidFluidOperator{2},(k,i,p)::Tuple{Any,Any,Any},(l,j,q)::Tuple{DataType,Any,Any}) = 
    _getindex_from_tuple_sign(B,(k,i,p),(l,j,q))
Base.getindex(B::FluidFluidOperator{2},(k,i,p)::Tuple{DataType,Any,Any},(l,j,q)::Tuple{DataType,Any,Any}) = 
    _getindex_from_tuple_sign(B,(k,i,p),(l,j,q))
# function Base.getindex(B::FluidFluidOperator{2},(k,i,p)::Tuple{Any,Any,Any},(l,j,q)::Tuple{Any,Any,Any})
#     k_map = findall(_sign_map(B,k))
#     l_map = findall(_sign_map(B,l))
#     (i==:)&&(i=1:n_phases(B.dq))
#     (j==:)&&(j=1:n_phases(B.dq))
#     (p==:)&&(p=1:n_bases_per_cell(B.dq))
#     (q==:)&&(q=1:n_bases_per_cell(B.dq))
#     k_map = k_map[map(x->x[1]∈i,k_map)]
#     l_map = l_map[map(x->x[1]∈j,l_map)]
#     return vcat([hcat([B[(k_map[row][2]-1,k_map[row][1],p),(l_map[col][2]-1,l_map[col][1],q)] for col in 1:length(l_map)]...) for row in 1:length(k_map)]...)
# end

Base.getindex(B::FluidFluidOperator{2},(k,i)::Tuple{Any,Any},(l,j)::Tuple{Any,Any}) = B[(k,i,:),(l,j,:)]
Base.getindex(B::FluidFluidOperator{2},k::Type{<:AbstractSign},l) = B[(k,:,:),(l,:,:)]
Base.getindex(B::FluidFluidOperator{2},k::Colon,l::Type{<:AbstractSign}) = B[(k,:,:),(l,:,:)]

index(ffq,s::Type{<:AbstractSign}) = _build_index_from_tuple_sign(ffq,(s,:,:))
