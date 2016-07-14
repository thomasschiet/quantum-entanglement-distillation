using Convex
using SCS

facts("Partial trace") do

  context("Function") do
    X_ab = reshape(collect(1:16), (4, 4))'
    X_a = [(1+6) (3+8); (9+14) (11+16)]
    X_b = [(1+11) (2+12); (5+15) (6+16)]

    @fact partialtrace(X_ab, 2, [2; 2]) --> X_a
    @fact partialtrace(X_ab, 1, [2; 2]) --> X_b

    X = rand(32, 32)
    @fact partialtrace(X, [1; 2; 3], [4; 4; 2])[1, 1] --> roughly(trace(X))
    @fact partialtrace(X, [2; 1; 3], [4; 4; 2])[1, 1] --> roughly(trace(X))
  end

  context("Convex.jl compatibility") do
    W = Semidefinite(4)
    problem = maximize(trace(partialtrace(W, 1, [2; 2])))
    problem.constraints += ([trace(W) ≤ 1])
    problem.constraints += ([partialtrace(W, 1, [2; 2]) == partialtrace(W, 2, [2; 2])])

    solve!(problem, SCSSolver(verbose = false))
    @fact vexity(problem) --> AffineVexity()
    @fact problem.optval --> roughly(1, TOL)
    @fact partialtrace(W.value, 1, [2; 2]) --> roughly(partialtrace(W.value, 2, [2; 2]), TOL)
  end
end
