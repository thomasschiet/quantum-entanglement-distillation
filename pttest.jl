# include("partialTrace.jl")
include("copies.jl")
include("eVec.jl")
include("wernerState.jl")
# include("quantumim.jl")

using Convex
using SCS

F = Semidefinite(256*2)
problem = minimize(trace(partialtrace(F, 2, [16, 16*2])))
# problem = minimize(trace(F))

problem.constraints += ([(eye(512) - F) in :SDP]);
solve!(problem,SCSSolver(verbose=true))
problem.optval
