# using FluidFluidQueues

##
# import Pkg
# Pkg.activate(".")
using DiscretisedFluidQueues, SparseArrays, LinearAlgebra, StableRNGs

# include("../src/operators.jl")
# include("../src/partition.jl")
##

using Test


T = [-2.5 2 0.5; 1 -2 1; 1 2 -3]
C = [0.0; 2.0; -6.0]

model = BoundedFluidQueue(T,C,6.0)
am = augment_model(model)
nodes = 0.0:1.0:6.0
order = 5
mesh = DGMesh(nodes,order)
dq = DiscretisedFluidQueue(am,mesh)
dq_notam = DiscretisedFluidQueue(model,mesh)
rates_lwr = [1.0; -3.0]
rates_upr = -[1.0; 2.0]
ffq_rates = Matrix(transpose([
    -1.0 0.0 3.0
    -1.0 2.0 0.0
    0.0 -2.0 3.0
    0.0 2.0 -3.0 
    1.0 -2.0 0.0
    1.0 0.0 -3.0
]))
# ffq_rates = [ones(1,length(nodes)-1);-2*ones(1,length(nodes)-1);-3*ones(1,length(nodes)-1)]
_rates = Rates(ffq_rates,rates_lwr,rates_upr) 
ffq = augment_model(FluidFluidQueue(dq_notam,_rates))
aug_rates = ffq.rates.interior
B = build_full_generator(dq)
ffq_B = FluidFluidGenerator(ffq)

for i in 1:n_phases(model)
    for j in 1:n_phases(model)
        @test B[(:,i,:),(:,j,:)]==ffq_B[(:,i),(:,j)]
        @test B[(:,i,:),(:,j,:)]==ffq_B[(:,i,:),(:,j,:)]
    end
end
R = RatesOperator(ffq)
for i in 1:n_phases(am)
    @show i
    lwr_idx = sum(DiscretisedFluidQueues._has_left_boundary(am.S)[1:i])
    upr_idx = sum(DiscretisedFluidQueues._has_right_boundary(am.S)[1:i])
    n₋ = DiscretisedFluidQueues._has_left_boundary(am.S,i) ? Int(rates_lwr[lwr_idx]==0.0) : 0
    n₊ = DiscretisedFluidQueues._has_right_boundary(am.S,i) ? Int(rates_upr[upr_idx]==0.0) : 0
    @test R[(Zero,i)]==zeros(sum(aug_rates[i,:].==0.0)*n_bases_per_cell(mesh)+n₊+n₋)
    n₋ = DiscretisedFluidQueues._has_left_boundary(am.S,i) ? Int(rates_lwr[lwr_idx]>0.0) : 0
    n₊ = DiscretisedFluidQueues._has_right_boundary(am.S,i) ? Int(rates_upr[upr_idx]>0.0) : 0
    @test R[(Plus,i)]==(1/max(1,i-1))*ones(sum(aug_rates[i,:].>0.0)*n_bases_per_cell(mesh)+n₊+n₋)
    n₋ = DiscretisedFluidQueues._has_left_boundary(am.S,i) ? Int(rates_lwr[lwr_idx]<0.0) : 0
    n₊ = DiscretisedFluidQueues._has_right_boundary(am.S,i) ? Int(rates_upr[upr_idx]<0.0) : 0
    @test R[(Minus,i)]==(1/max(1,i-1))*ones(sum(aug_rates[i,:].<0.0)*n_bases_per_cell(mesh)+n₊+n₋)
end
p_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(Plus,:,:)))
m_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(Minus,:,:)))
z_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(Zero,:,:)))
np_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotPlus,:,:)))
nm_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotMinus,:,:)))
nz_idx = Set(FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotZero,:,:)))
@test intersect(p_idx,m_idx)==Set{Int64}()
@test intersect(z_idx,m_idx)==Set{Int64}()
@test intersect(z_idx,p_idx)==Set{Int64}()
@test intersect(nz_idx,z_idx)==Set{Int64}()
@test intersect(np_idx,p_idx)==Set{Int64}()
@test intersect(nm_idx,m_idx)==Set{Int64}()
@test intersect(nm_idx,p_idx)==p_idx
@test intersect(nm_idx,z_idx)==z_idx
@test intersect(np_idx,m_idx)==m_idx
@test intersect(np_idx,z_idx)==z_idx
@test intersect(nz_idx,p_idx)==p_idx
@test intersect(nz_idx,m_idx)==m_idx
@test union(p_idx,m_idx,z_idx)==Set(1:length(R))

@test R[Zero]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(Zero,:,:))]
@test R[Plus]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(Plus,:,:))]
@test R[Minus]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(Minus,:,:))]
@test R[NotZero]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotZero,:,:))]
@test R[NotPlus]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotPlus,:,:))]
@test R[NotMinus]==R[FluidFluidQueues._build_index_from_tuple_sign(ffq,(NotMinus,:,:))]

D = InOutGenerator(ffq_B,0.0)
psi = build_psi(D)
xi = build_xi(ffq_B,psi)
marginal, point_mass, K = build_limit_dist_operators(ffq_B,D,R,psi,xi)

macro resetrng1()
    return :(StableRNG(04012022))
end
macro resetrng2()
    return :(StableRNG(204012022))
end
rng1 = @resetrng1()
n_sim = 100
XYvals = 6*rand(rng1,n_sim)
phi0 = rand(rng1,1:4,n_sim)
interior_rates = repeat(C,1,length(nodes)-1)
ffq = augment_model(FluidFluidQueue(dq_notam,Rates(interior_rates,[0.0;0.0],[0.0;0.0])))
ffq_notam = FluidFluidQueue(dq_notam,Rates(interior_rates,[0.0;0.0],[0.0;0.0]))
rng2 = @resetrng2()
SFM, YSims = simulate(ffq,FluidFluidQueues.n_jumps_y(10),
    (φ=phi0, X=XYvals, Y=XYvals),
    rng2,
)
rng2 = @resetrng2()
SFM_notam, YSims_notam = simulate(ffq_notam,FluidFluidQueues.n_jumps_y(10),
    (φ=max.(1,phi0.-1), X=XYvals, Y=XYvals), 
    rng2,
)

@test SFM.X ≈ YSims
@test SFM_notam.X ≈ YSims_notam
@test SFM_notam.t ≈ SFM.t
@test SFM_notam.φ ≈ max.(1,SFM.φ.-1)