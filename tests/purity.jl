facts("Purity") do

  W = eye(4)
  t = 1

  function findW(X_opt, W_opt, n::Int = 1)
    N = size(X)[1]
    W = Semidefinite(N)
    w = trace(X_opt * W_opt)
    # w = 10
    # println(w)

    problem = minimize(trace(X_opt) + w * trace(X_opt*W))
    problem.constraints += 0 ⪯ W
    problem.constraints += W ⪯ eye(N)
    problem.constraints += trace(W) == N - n

    solve!(problem, SCSSolver(verbose = false))
    return (W.value, problem.optval)
  end

  function findX(X_opt, W_opt, t)
    X = Semidefinite(4)
    w = trace(X_opt * W_opt)
    # w = 10
    # println(w)

    problem = minimize(trace(X) + w * trace(W_opt*X))
    problem.constraints += [
        trace(X) ≥ 1
        X[2, 1] == 1/2
        # trace(X * W_opt) ≤ t
    ]
    solve!(problem, SCSSolver(verbose = false))
    return X.value
  end

  W = rand(4, 4)
  W = W - W'
  X = eye(4)
  t = 1
  i = 1
  function ispure(X::AbstractMatrix; tol::Number = 0)
    if tol == 0
      return rank(X) == 1
    else
      return ispure(eigvals(X))
    end
  end

  function ispure(λs::Vector; tol::Number = 0)
    λs = abs(λs)
    λs = -sort(-λs)
    return λs[2] ≤ tol && λs[1] ≥ tol
  end

  while i < 1000
    X = findX(X, W, t)

    if ispure(X, 1e-6)
      break
    end

    (W, t) = findW(X, W, 1)
  end

  X = round(X, 3)

  @fact rank(X) --> 1
end
