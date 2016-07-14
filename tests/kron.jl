using Convex
using SCS

facts("Kron") do
  context("Convex.jl compatibility") do
    X = Semidefinite(4)
    id = Constant(eye(4))

    W = id ⊗ X

    problem = maximize(trace(W))
    problem.constraints += trace(X) ≤ 1

    solve!(problem, SCSSolver(verbose = false))

    @fact problem.optval --> roughly(4, TOL)
  end
end
