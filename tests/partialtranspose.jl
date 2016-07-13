using Convex
using SCS

facts("Partial transpose") do
  context("Naive") do
    prePT = reshape(1:16, 4, 4)'
    postPT = [1 5 3 7; 2 6 4 8; 9 13 11 15; 10 14 12 16]
    @fact partialtranspose(prePT, 2, naive = true) --> postPT
  end

  context("QETLAB style") do
    prePT = reshape(1:16, 4, 4)'
    postPT = [1 5 3 7; 2 6 4 8; 9 13 11 15; 10 14 12 16]
    @fact partialtranspose(prePT, 2, [2; 2]) --> postPT

    prePT = reshape(1:16, 4, 4)'
    postPT = [1 2 9 10; 5 6 13 14; 3 4 11 12; 7 8 15 16];
    @fact partialtranspose(prePT, 1, [2; 2]) --> postPT

    prePT = rand(4, 4)
    @fact partialtranspose(partialtranspose(prePT, 1, [2; 2]), 2, [2; 2]) --> prePT'
    @fact partialtranspose(prePT, [1; 2], [2; 2]) --> prePT'
  end

  context("Convex.jl compatibility") do
    W = Semidefinite(4)
    id = eye(4)

    problem = maximize(trace(W))
    problem.constraints = [
      W ⪯ id
      W[2, 2] == 1
      W == partialtranspose(W, 2)
    ]

    solve!(problem, SCSSolver(verbose = false))

    @fact problem.optval --> roughly(4, TOL)
    @fact partialtranspose(W.value, 2) --> roughly(W.value, TOL)
  end
end
