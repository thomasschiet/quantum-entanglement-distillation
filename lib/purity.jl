using Convex
using SCS

export purity

"""
For some input state ρ, find an optimal Choi state that maximizes fidelity with rank 1. (purity condition)
"""
function purity(ρ::AbstractArray)
  nAhat = 2
  nF_A = 2
  nAprime = 4
  nE_A = 2
  nBhat = nAhat
  nF_B = nF_A
  nBprime = nAprime
  nE_B = nE_A

  nA = nAhat * nF_A * nAprime * nE_A
  nB = nBhat * nF_B * nBprime * nE_B
  nC = nA * nB

  t_A = 1
  t_B = 1

  # instantiate direction matrices
  W_Abar = eye(nA)
  W_Bbar = eye(nB)

  # TODO: loop this properly
  for i = 1:10
    (C_AbarBbar, ) = findChoi(ρ, W_Abar, W_Bbar, t_A, t_B)
    (W_Abar, W_Bbar, t_A, t_B) = findWs(C_AbarBbar, )
  end

  return C_AbarBbar
end

"""
Find a Choi state given direction matrices
"""
function findChoi(ρ::AbstractArray, W_Abar::AbstractArray, W_Bbar::AbstractArray, t_A::Number, t_B::Number; verbose::Bool = false)
  @assert isQuantumState(ρ) "expected ρ to be a quantum state"
  @assert issquare(W_Abar) "expected W_Abar to be square"
  @assert issquare(W_Bbar) "expected W_Bbar to be square"
  @assert 0 ≤ t_A ≤ 1
  @assert 0 ≤ t_B ≤ 1

  nAhat = 2
  nF_A = 2
  nAprime = 4
  nE_A = 2
  nBhat = nAhat
  nF_B = nF_A
  nBprime = nAprime
  nE_B = nE_A

  nA = nAhat * nF_A * nAprime * nE_A
  nB = nBhat * nF_B * nBprime * nE_B
  nC = nA * nB

  C_AbarBbar = Semidefinite(nC)
  C_Abar = partialtrace(C_AB, 2, [nA, nB])
  C_Bbar = partialtrace(C_AB, 1, [nA, nB])

  problem = maximize(nA*nB*trace(ρ' * C_AbarBbar))

  # constraints from purity
  problem.constraints += ([(trace(C_Abar * W_Abar)) ≤ t_A])
  problem.constraints += ([(trace(C_Bbar * W_Bbar)) ≤ t_B])

  solve!(problem, SCSSolver(verbose = verbose))

  return (C_AbarBbart.value, )
end

"""
Find a direction matrix given a Choi state for a specific system
"""
function findW(C_opt::AbstractArray, system::Union{AbstractString, Int}; verbose::Bool = false)
  @assert system == "A" || system == "B" || system == 1 || system == 2 "expected system to be A or B"
  @assert isQuantumState(C_opt) "expected C_opt to be a quantum state"

  if system == "A"
    system = 1
  elseif system == "B"
    system = 2
  end

  nC = size(C_opt)[1]
  nX = round(Int, sqrt(nC))
  C_X = partialtrace(C_opt, system, [nX; nX])

  W_Xbar = Variable(size(C_X))

  problem = minimize(trace(C_X' * W_Xbar))
  problem.constraints += ([0 ⪯ W_Xbar ⪯ eye(size(W_Xbar))])
  problem.constrains += ([trace(W_Xbar) = size(C_X)[1] - 1])

  solve!(SCSSolver(verbose = verbose))
  return problem.optval
end

"""
Find W for A and for B
"""
function findWs(C_opt::AbstractArray; verbose::Bool = false)
  @assert issquare(C_opt) "expected C_opt to be square"

  (W_Abar, t_A) = findW(C_opt, "A", verbose = verbose)
  (W_Bbar, t_B) = findW(C_opt, "B", verbose = verbose)
  return (W_Abar, W_Bbar, t_A, t_B)
end
