export constructNet
export constructRadialNet
export findChoiStateA
export findChoiStateB

using Iterators
using Convex
using SCS

function constructNet(dim::Int, steps::Int)
  N = dim # dimensions
  m = steps # steps

  region = -1:(2/(steps-1)):1

  return map((x) -> collect(x), collect(product(ntuple(x -> region, dim)...)))
end

function constructRadialNet(dim::Int, steps::Int, factor::Number)
  N = dim # dimensions
  m = steps # steps

  # exclude the zero
  r = collect(linspace(0, 1, steps+1))[2:end]
  r = [1]

  # exclude 2π
  Θ = collect(linspace(0, 2, 2*steps+1))[1:end-1]
  Θ *= π

  region = collect(linspace(0, 1, steps+2))[2:end-1]
  region = collect(linspace(0, 1, steps))
  region *= π

  if dim-2 > 0
    parameters = map((x) -> collect(x), collect(product(r, Θ, ntuple(x -> region, dim-2)...)))
  else
    parameters = map((x) -> collect(x), collect(product(r, Θ)))
  end

  # map the spherical coordinates to cartesian coordinates
  return map((p) -> sphericalToCartesian(p, dim)*factor, parameters)
end

function sphericalToCartesian(p::Vector, dim::Int)
  x = p
  # set all coordinates of x
  for i in 1:dim
    # set the coordinate equal to r
    x[i] = p[1]

    # now multiply with the trig functions
    # i=1, j=2..2 x[1] = r cos(p[2])
    # i=2  j=2..3 x[2] = r sin(p[2])cos(p[3])
    for j in 2:(i+1)
      if j > dim
        continue
      end

      # the last factor should be a cosine
      if j == i+1
        x[i] *= cos(p[j])
      else
        x[i] *= sin(p[j])
      end
    end
  end

  return x
end

#
function findChoiState(c::Vector, basis::Array{Any, 1}, ρ::Union{Matrix, SparseMatrixCSC}, δ::Number, ϵ::Number, state::AbstractString = "A"; verbose::Bool = false)
  ps = Convex.Variable[]
  qs = Convex.Variable[]

  # define the objective
  # the p and q variables are needed to make the problem convex
  # this works because if we want to minimize |x - a| st x ∈ X
  # this can be converted to
  #
  # min  p + q
  # s.t. p - q = x - a
  #      p, q ≥ 0
  #      x ∈ X
  objective = 0
  for b in basis
    append!(ps, Convex.Variable[Variable()])
    append!(qs, Convex.Variable[Variable()])
    objective += ps[end] + qs[end]
  end
  problem = minimize(objective)

  # define Choi states
  C_size = round(Int, 2 * sqrt(size(ρ)[1]))
  C_BhatBprime = Semidefinite(C_size)
  C_0 = Semidefinite(C_size)
  C_1 = Semidefinite(C_size)

  # add constraints for absolute values
  i = 1
  for b in basis
    # problem.constraints += ([ps[i] - qs[i] == c[i] - trace(ptrace(b, 1, [2; 3]) * ptrace(C_BhatBprime, 1, [2; 3]))])
    problem.constraints += ([ps[i] - qs[i] == c[i] - trace(b * C_BhatBprime)])
    problem.constraints += ([ps[i] ≥ 0])
    problem.constraints += ([qs[i] ≥ 0])
    i += 1
  end

  # add probability constraints dependent on A/B
  if state == "A"
    problem.constraints += ([3 * trace(ptrace(ρ, 2, [3; 3]) * ptrace(C_1, 1, [2; 3])) ≥ δ - ϵ])
    problem.constraints += ([3 * trace(ptrace(ρ, 2, [3; 3]) * ptrace(C_1, 1, [2; 3])) ≤ δ + ϵ])
  else
    problem.constraints += ([3 * trace(ptrace(ρ, 1, [3; 3]) * ptrace(C_1, 1, [2; 3])) ≥ δ - ϵ])
    problem.constraints += ([3 * trace(ptrace(ρ, 1, [3; 3]) * ptrace(C_1, 1, [2; 3])) ≤ δ + ϵ])
  end

  # the following constraints are to ensure C is a Choi state
  problem.constraints += ([C_BhatBprime == C_1 + C_0])
  problem.constraints += ([ptrace(C_BhatBprime, 1, [2; 3]) == eye(3)/3])
  problem.constraints += ([C_0 ⪯ eye(C_size)/C_size])
  problem.constraints += ([C_1 ⪯ eye(C_size)/C_size])
  solve!(problem, SCSSolver(verbose = verbose))
  if verbose
    if abs(problem.optval) ≤ 10e-6
      println(string("🎉 found Choi state at a distance of ", problem.optval, "≈ 0"), :green)
    else
      println(string("❌ found Choi state at a distance of ", problem.optval), :red)
    end
    println(string(round(C_BhatBprime.value, 3)))

    i = 1
    for p in ps
      x = p.value + qs[i].value
      if x > 10e-6
        println(string("i = " , i, " equals", x), [:bold, :red, :bgWhite])
      end
      i += 1
    end
  end
  return (problem, C_BhatBprime, C_1, C_0)
end

# aliasses
constructRadialNet(dim::Int, steps::Int) = constructRadialNet(dim, steps, 1)
findChoiStateA(c, basis::Array{Any, 1}, ρ::Union{Matrix, SparseMatrixCSC}, δ::Number, ϵ::Number) = findChoiState(c, basis, ρ, δ, ϵ, "A")
findChoiStateA(c, basis::Array{Any, 1}, ρ::Union{Matrix, SparseMatrixCSC}, δ::Number, ϵ::Number; verbose::Bool = false) = findChoiState(c, basis, ρ, δ, ϵ, "A", verbose = verbose)
findChoiStateB(c, basis::Array{Any, 1}, ρ::Union{Matrix, SparseMatrixCSC}, δ::Number, ϵ::Number) = findChoiState(c, basis, ρ, δ, ϵ, "B")
findChoiStateB(c, basis::Array{Any, 1}, ρ::Union{Matrix, SparseMatrixCSC}, δ::Number, ϵ::Number; verbose::Bool = false) = findChoiState(c, basis, ρ, δ, ϵ, "B", verbose = verbose)
