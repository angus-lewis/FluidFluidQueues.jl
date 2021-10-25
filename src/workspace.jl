import Pkg
Pkg.activate(".")
using DiscretisedFluidQueues

include("src/operators.jl")
include("src/partition.jl")

T = [-2.5 2 0.5; 1 -2 1; 1 2 -3]
C = [0.0; 2.0; -3.0]

model = BoundedFluidQueue(T,C)

nodes = 0.0:1.0:6.0
order = 2
mesh = DGMesh(nodes,order)
dq = DiscretisedFluidQueue(model,mesh)
B = build_lazy_generator(dq);
rates_lwr = [1.0; -1.0]
rates_upr = -[1.0; 1.0]
ffq_rates = Matrix(transpose([
    -1.0 0.0 1.0
    -1.0 1.0 0.0
    0.0 -1.0 1.0
    0.0 1.0 -1.0 
    1.0 -1.0 0.0
    1.0 0.0 -1.0
]))
ffq = FluidFluidQueue(B,ffq_rates,rates_lwr,rates_upr);
R = RatesOperator(ffq)
R[Plus]
ffq[(:,2), (:,2)]