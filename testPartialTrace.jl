include("BEPQuantum.jl")

using Convex
using SCS
using .BEPQuantum

workspace()
dims = [4; 4; 4; 4; 4; 3]
println(prod(dims))
W = Semidefinite(prod(dims))

problem = minimize(trace(W))
problem.constraints += ([W[3, 1] == 1])
problem.constraints += ([W[4, 1] == 1])
problem.constraints += ([partialtrace(W, 1, dims) == partialtrace(W, 2, dims)])
problem.constraints += ([partialtrace(W, 3, dims) == partialtrace(W, 4, dims)])
# problem.constraints += ([partialtrace(W, 5, dims) == partialtrace(W, 6, dims)])
@time solve!(problem, SCSSolver(verbose = true))
println(round(problem.optval, 2))
println(round(W.value, 2))
# partialtrace(W, 1, [4; 4])
println(round(partialtrace(W.value, 1, [4; 4]), 2))
println(round(partialtrace(W.value, 2, [4; 4]), 2))


traceOut = 3 
traceOut_i = 1
dims = [traceOut; 2]
totalDim = prod(dims)
println(totalDim)
X = reshape(1:totalDim^2, (totalDim, totalDim))'
map = zeros(totalDim^2, round(Int, (totalDim / traceOut)^2))
for i = 1:totalDim^2
  E = zeros(totalDim^2)
  E[i] = 1
  E = reshape(E, (totalDim, totalDim))
  map[i, :] = vec(partialtrace(E, traceOut_i, dims))
end

map = sparse(map)
println(X)
println(map)

@show map






p

X
reshape(map'*vec(X), (2, 2))
