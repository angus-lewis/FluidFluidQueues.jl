###############
###############
## partition ##
###############
###############

""" First index of the generator to index the generator with ("+","-") etc, elements are strings """
const PlusMinusIndex = Union{String,Tuple{String,String}}

""" Second index of the generator to index the generator with (i,j) etc, elements are Int or :"""
const PhaseIndex = Union{Tuple{Union{Int64,Colon},Union{Int64,Colon}},Int64,Colon}

""" Index of the generator e.g. ("+","-"),(i,j) etc, elements are (String,String),(Int/Colon,Int/Colon) """
const GeneratorIndex = Union{PlusMinusIndex,Tuple{PlusMinusIndex,PhaseIndex}}

""" The dictionary Fil which indicate which cells correspond to +, -, or 0 fluid-fluid rates """
const IndexDict = Dict{Tuple{String,Union{Int64,Colon}},BitArray{1}}

""" The dictionary which returns various partitions of the generator i.e. ("+","-"),(i,j) """
const PartitionedGenerator = Dict{GeneratorIndex,SparseArrays.SparseMatrixCSC{Float64,Int64}}

"""

    MakeFil(
        model::FluidFluidQueue,
        Nodes::Array{<:Real,1},
        )

Construct dict with entries indexing which cells belong to Fᵢᵐ. 
"""
function MakeFil(
    model::FluidFluidQueue,
    Nodes::Array{<:Real,1},
    )
    meshNIntervals = length(Nodes) - 1
    Δtemp = Nodes[2:end] - Nodes[1:end-1]

    Fil = IndexDict()
    
    ## Construct the sets Fᵐ = ⋃ᵢ Fᵢᵐ, global index for sets of type m
    idxPlus = model.r.r(Nodes[1:end-1].+Δtemp[:]/2).>0
    idxZero = model.r.r(Nodes[1:end-1].+Δtemp[:]/2).==0
    idxMinus = model.r.r(Nodes[1:end-1].+Δtemp[:]/2).<0
    for i in 1:NPhases(model)
        Fil[("+",i)] = idxPlus[:,i]
        Fil[("0",i)] = idxZero[:,i]
        Fil[("-",i)] = idxMinus[:,i]
        if model.C[i] .<= 0
            Fil[("p+",i)] = [model.r.r(model.Bounds[1,1])[i]].>0
            Fil[("p0",i)] = [model.r.r(model.Bounds[1,1])[i]].==0
            Fil[("p-",i)] = [model.r.r(model.Bounds[1,1])[i]].<0
        end
        if model.C[i] .>= 0
            Fil[("q+",i)] = [model.r.r(model.Bounds[1,end])[i]].>0
            Fil[("q0",i)] = [model.r.r(model.Bounds[1,end])[i]].==0
            Fil[("q-",i)] = [model.r.r(model.Bounds[1,end])[i]].<0
        end
    end
    currKeys = keys(Fil)
    for ℓ in ["+", "-", "0"], i = 1:NPhases(model)
        if ((ℓ,i) ∉ currKeys)
            Fil[(ℓ,i)] = falses(meshNIntervals)
        end
        if (("p"*ℓ,i) ∉ currKeys) && (model.C[i] <= 0)
            Fil["p"*ℓ,i] = falses(1)
        end
        if (("p"*ℓ,i) ∉ currKeys) && (model.C[i] > 0)
            Fil["p"*ℓ,i] = falses(0)
        end
        if (("q"*ℓ,i) ∉ currKeys) && (model.C[i] >= 0)
            Fil["q"*ℓ,i] = falses(1)
        end
        if (("q"*ℓ,i) ∉ currKeys) && (model.C[i] < 0)
            Fil["q"*ℓ,i] = falses(0)
        end
    end
    for ℓ in ["+", "-", "0"]
        Fil[(ℓ,:)] = falses(meshNIntervals * NPhases(model))
        Fil[("p"*ℓ,:)] = trues(0)
        Fil[("q"*ℓ,:)] = trues(0)
        for i = 1:NPhases(model)
            idx = findall(Fil[(ℓ,i)]) .+ (i - 1) * meshNIntervals
            Fil[(ℓ,:)][idx] .= true
            Fil[("p"*ℓ,:)] = [Fil[("p"*ℓ,:)]; Fil[("p"*ℓ,i)]]
            Fil[("q"*ℓ,:)] = [Fil[("q"*ℓ,:)]; Fil[("q"*ℓ,i)]]
        end
    end
    return Fil
end

function MakeDict(
    B::Union{Array{<:Real,2},SparseArrays.SparseMatrixCSC{<:Real,Int64}},
    model::Model, 
    mesh::Mesh,
    Fil::IndexDict,
    )

    ## Make a Dictionary so that the blocks of B are easy to access
    # N₋ = sum(model.C.<=0)
    # N₊ = sum(model.C.>=0)

    BDict = PartitionedGenerator()

    ppositions = cumsum(model.C .<= 0)
    qpositions = cumsum(model.C .>= 0)
    for ℓ in ["+", "-", "0"], m in ["+", "-", "0"]
        for i = 1:NPhases(model), j = 1:NPhases(model)
            FilBases = repeat(Fil[(ℓ,i)]', NBases(mesh), 1)[:]
            pitemp = falses(N₋(model.S))
            qitemp = falses(N₊(model.S))
            pjtemp = falses(N₋(model.S))
            qjtemp = falses(N₊(model.S))
            if model.C[i] <= 0
                if length(pitemp) > 0 
                    pitemp[ppositions[i]] = Fil[("p"*ℓ,i)][1]
                end
            end
            if model.C[j] <= 0
                if length(pjtemp) > 0
                    pjtemp[ppositions[j]] = Fil[("p"*m,j)][1]
                end
            end
            if model.C[i] >= 0
                if length(qitemp) > 0
                    qitemp[qpositions[i]] = Fil[("q"*ℓ,i)][1]
                end
            end
            if model.C[j] >= 0
                if length(qjtemp) > 0
                    qjtemp[qpositions[j]] = Fil[("q"*m,j)][1]
                end
            end
            i_idx = [
                pitemp
                falses((i - 1) * TotalNBases(mesh))
                FilBases
                falses(NPhases(model) * TotalNBases(mesh) - i * TotalNBases(mesh))
                qitemp
            ]
            FjmBases = repeat(Fil[(m, j)]', NBases(mesh), 1)[:]
            j_idx = [
                pjtemp
                falses((j - 1) * TotalNBases(mesh))
                FjmBases
                falses(NPhases(model) * TotalNBases(mesh) - j * TotalNBases(mesh))
                qjtemp
            ]
            BDict[((ℓ, m),(i, j))] = B[i_idx, j_idx]
        end
        # below we need to use repeat(Fil[ℓ]', NBases(mesh), 1)[:] to
        # expand the index Fil[ℓ] from cells to all basis function
        FlBases =
            [Fil["p"*ℓ,:]; repeat(Fil[ℓ,:]', NBases(mesh), 1)[:]; Fil["q"*ℓ,:]]
        FmBases =
            [Fil["p"*m,:]; repeat(Fil[m,:]', NBases(mesh), 1)[:]; Fil["q"*m,:]]
        BDict[((ℓ,m),(:,:))] = B[FlBases, FmBases]
    end
    return BDict
end

function MakeDict(
    B::Union{Array{<:Real,2},SparseArrays.SparseMatrixCSC{<:Real,Int64}},
    C::Array{<:Real,1},
    order::Int,
    Fil::IndexDict,
    )

    ## Make a Dictionary so that the blocks of B are easy to access
    N₋ = sum(C.<=0)
    N₊ = sum(C.<=0)

    BDict = PartitionedGenerator()

    n_phases = length(C)
    total_n_bases = Int((size(B,1)-N₋-N₊)/n_phases)    

    ppositions = cumsum(C .<= 0)
    qpositions = cumsum(C .>= 0)
    for ℓ in ["+", "-", "0"], m in ["+", "-", "0"]
        for i = 1:n_phases, j = 1:n_phases
            FilBases = repeat(Fil[(ℓ,i)]', order, 1)[:]
            pitemp = falses(N₋)
            qitemp = falses(N₊)
            pjtemp = falses(N₋)
            qjtemp = falses(N₊)
            if C[i] <= 0
                if length(pitemp) > 0 
                    pitemp[ppositions[i]] = Fil[("p"*ℓ,i)][1]
                end
            end
            if C[j] <= 0
                if length(pjtemp) > 0
                    pjtemp[ppositions[j]] = Fil[("p"*m,j)][1]
                end
            end
            if C[i] >= 0
                if length(qitemp) > 0
                    qitemp[qpositions[i]] = Fil[("q"*ℓ,i)][1]
                end
            end
            if C[j] >= 0
                if length(qjtemp) > 0
                    qjtemp[qpositions[j]] = Fil[("q"*m,j)][1]
                end
            end
            i_idx = [
                pitemp
                falses((i - 1) * total_n_bases)
                FilBases
                falses(n_phases * total_n_bases - i * total_n_bases)
                qitemp
            ]
            FjmBases = repeat(Fil[(m, j)]', order, 1)[:]
            j_idx = [
                pjtemp
                falses((j - 1) * order)
                FjmBases
                falses(n_phases * order - j * order)
                qjtemp
            ]
            BDict[((ℓ, m),(i, j))] = B[i_idx, j_idx]
        end
        # below we need to use repeat(Fil[ℓ]', order, 1)[:] to
        # expand the index Fil[ℓ] from cells to all basis function
        FlBases =
            [Fil["p"*ℓ,:]; repeat(Fil[ℓ,:]', order, 1)[:]; Fil["q"*ℓ,:]]
        FmBases =
            [Fil["p"*m,:]; repeat(Fil[m,:]', order, 1)[:]; Fil["q"*m,:]]
        BDict[((ℓ,m),(:,:))] = B[FlBases, FmBases]
    end
    return BDict
end

##############
##############
## indexing ##
##############
##############

# function getindex_correct(B::LazyGenerator,row::Int,col::Int)
#     checkbounds(B,row,col)

#     ei = zeros(1,size(B,1))
#     ei[row] = 1
#     return (ei*B)[col]
# end

function getindex(B::LazyGenerator,plus_minus_index::PlusMinusIndex,phase_index::PhaseIndex)
    N₋ = sum(B.C.<=0)
    N₊ = sum(B.C.<=0)

    out = []

    ppositions = cumsum(B.C .<= 0)
    qpositions = cumsum(B.C .>= 0)

    ℓ,m = plus_minus_index
    if typeof(phase_index[1])==Colon
        i_range = 1:length(B.C)
    else
        i_range = phase_index[1]
    end
    if typeof(phase_index[2])==Colon
        j_range = 1:length(B.C)
    else
        j_range = phase_index[2]
    end

    for i = i_range, j = j_range
        FilBases = repeat(B.Fil[ℓ,i]', size(B.D,1), 1)[:]
        pitemp = falses(N₋)
        qitemp = falses(N₊)
        pjtemp = falses(N₋)
        qjtemp = falses(N₊)
        if B.C[i] <= 0
            if length(pitemp) > 0 
                pitemp[ppositions[i]] = B.Fil["p"*ℓ,i][1]
            end
        end
        if B.C[j] <= 0
            if length(pjtemp) > 0
                pjtemp[ppositions[j]] = B.Fil["p"*m,j][1]
            end
        end
        if B.C[i] >= 0
            if length(qitemp) > 0
                qitemp[qpositions[i]] = B.Fil["q"*ℓ,i][1]
            end
        end
        if B.C[j] >= 0
            if length(qjtemp) > 0
                qjtemp[qpositions[j]] = B.Fil["q"*m,j][1]
            end
        end
        i_idx = [
            pitemp
            falses((i - 1) * size(B.D,1)*length(B.Δ))
            FilBases
            falses(length(B.C) * size(B.D,1)*length(B.Δ) - i * size(B.D,1)*length(B.Δ))
            qitemp
        ]
        FjmBases = repeat(B.Fil[m,j]', size(B.D,1), 1)[:]
        j_idx = [
            pjtemp
            falses((j - 1) * size(B.D,1)*length(B.Δ))
            FjmBases
            falses(length(B.C) * size(B.D,1)*length(B.Δ) - j * size(B.D,1)*length(B.Δ))
            qjtemp
        ]
        if (typeof(phase_index[1]) != Colon) && (typeof(phase_index[2]) != Colon)
            out = B[i_idx, j_idx]
        end
    end
    # we use repeat(mesh.Fil[ℓ]', NBases(mesh), 1)[:] to
    # expand the index mesh.Fil[ℓ] from cells to all basis function
    if (typeof(phase_index[1]) == Colon) && (typeof(phase_index[2]) == Colon)
        FlBases =
            [B.Fil["p"*ℓ,:]; repeat(B.Fil[ℓ,:]', size(B.D,1), 1)[:]; B.Fil["q"*ℓ,:]]
        FmBases =
            [B.Fil["p"*m,:]; repeat(B.Fil[m,:]', size(B.D,1), 1)[:]; B.Fil["q"*m,:]]
        out = B[FlBases, FmBases]
    end

    return out
end

getindex(B::LazyGenerator,idx::Tuple{PlusMinusIndex,PhaseIndex}) = getindex(B,idx[1],idx[2])

# this is a bit slow...
function MakeDict(B::LazyGenerator)
    BDict = PartitionedGenerator()
    for i in 1:length(B.C), j in 1:length(B.C)
        for k in ["+","-","0"], ℓ in ["+","-","0"]
            plus_minus_index = (k,ℓ)
            phase_index = (i,j)
            BDict[plus_minus_index, phase_index] = B[plus_minus_index, phase_index]
        end
    end
    i = :;
    j = :;
    for k in ["+","-","0"], ℓ in ["+","-","0"] 
        plus_minus_index = (k,ℓ)
        phase_index = (i,j)
        BDict[plus_minus_index, phase_index] = B[plus_minus_index, phase_index]
    end
    return BDict
end