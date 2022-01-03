import Pkg
Pkg.activate(".")
using DiscretisedFluidQueues
using LinearAlgebra
using Random 
using SparseArrays

include("queue.jl")
include("operators.jl")
include("partition.jl")
include("simulate.jl")

# include("FluidFluidQueues.jl")

# using FluidFluidQueues

T = [-2.5 2 0.5; 1 -2 1; 1 2 -3]
C = [0.0; 2.0; -3.0]

model = BoundedFluidQueue(T,C,6.0)

nodes = 0.0:1.0:6.0
order = 2
mesh = DGMesh(nodes,order)
dq = DiscretisedFluidQueue(model,mesh)
B = build_lazy_generator(dq);
rates_lwr = [0.0; 0.0]
rates_upr = [0.0; 0.0]
ffq_rates = Matrix(transpose([
    zeros(6) fill(2.0,6) fill(-3.0,6)
]))
ffq = FluidFluidQueue(B.dq,ffq_rates,rates_lwr,rates_upr);
R = RatesOperator(ffq)
R[Plus]

SFM, YSims = simulate(ffq,first_exit_y(eps(),6),(φ=ones(Int,100),X=fill(eps(),100),Y=fill(eps(),100)),MersenneTwister(1))