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

  context("Linear transformation") do
    function kronT(X::AbstractMatrix)
      e = AbstractMatrix[]
      push!(e, [1 0; 0 0])
      push!(e, [0 0; 1 0])
      push!(e, [0 1; 0 0])
      push!(e, [0 0; 0 1])

      T = zeros(16, 4)
      for i = 1:4
        T[:, i] = vec(kron(e[i], X))
      end
      return T
    end

    function kron2(x, y)
      T = kronT(y)
      return reshape(T * vec(x), (4, 4))
    end

    a = rand(2, 2)
    b = rand(2, 2)

    @fact kron2(a, b) --> kron(a, b)

    X = Semidefinite(2)
    id = eye(2)

    W = kron2(X, id)

    problem = maximize(trace(W))
    problem.constraints += trace(X) ≤ 1

    solve!(problem, SCSSolver(verbose = false))

    @fact problem.optval --> roughly(2, TOL)
  end
end
