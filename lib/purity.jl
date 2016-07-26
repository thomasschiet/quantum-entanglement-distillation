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
  W_Abar = zeros(nA, nA)
  W_Bbar = zeros(nB, nB)

  # do iterations until purity is reached or 100 iterations were done
  for i = 1:100
    println("Starting iteration ", i)
    @time (C_AbarBbar, problem, F, p_succ) = findChoi(ρ, W_Abar, W_Bbar, t_A, t_B)

    # sometimes Convex.jl fails.
    # try again in that case
    j = 1
    while C_AbarBbar == zeros((1024, 1024))
      println(" trying again ", j, "...")
      j += 1
      @time (C_AbarBbar, problem, F, p_succ) = findChoi(ρ, W_Abar, W_Bbar, t_A, t_B)
    end

    λ = eigvals(C_AbarBbar)
    λ = -sort(-real(λ))
    println(" λs = ", round(λ[1:10], 2))
    println(" F = ", F)
    println(" p_succ = ", p_succ)

    # check with some tolerance if the state is pure
    # if so, return the found Choi state
    if λ[2] ≤ 1e-3
      return C_AbarBbar
    end

    @time (W_Abar, W_Bbar, t_A, t_B) = findWs(C_AbarBbar, )
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

  δ = .5

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

  e0 = eVec(2, 1)
  e1 = eVec(2, 2)

  e11 = e1 ⊗ e1
  e11 = e11 * e11'

  epr = e0 ⊗ e0 + e1 ⊗ e1
  epr = epr * epr'

  e1 = e1 * e1'

  C_AbarBbar = Semidefinite(nC)
  C_Abar = partialtrace(C_AbarBbar, 2, [nA; nB])
  C_Bbar = partialtrace(C_AbarBbar, 1, [nA; nB])

  C_AprimeBprime = partialtrace(C_AbarBbar, [1; 2; 4; 5; 6; 8] , [nAhat; nF_A; nAprime; nE_A; nBhat; nF_B; nBprime; nE_B])

  problem = maximize(trace(epr ⊗ e11 ⊗ eye(nE_A * nE_B) ⊗ ρ' * C_AbarBbar))

  # constraints from purity
  problem.constraints += [trace(C_Abar * W_Abar) ≤ t_A]
  problem.constraints += [trace(C_Abar * W_Bbar) ≤ t_B]

  # constraints for Choi state
  problem.constraints += [C_AprimeBprime == eye(16)/16]
  problem.constraints += [trace(C_AbarBbar) == 1]

  # PPT constraint
  problem.constraints += [ptranspose(C_AbarBbar) == C_AbarBbar]

  # probability constraints
  problem.constraints += [trace(eye(nAhat) ⊗ e1 ⊗ partialtrace(ρ, 2, [4; 4]) ⊗ eye(nE_A) ⊗ eye(nBhat) ⊗ e1 ⊗ partialtrace(ρ, 1, [4; 4]) ⊗ eye(nE_B) * C_AbarBbar) ≤ δ]
  problem.constraints += [trace(eye(nAhat) ⊗ e1 ⊗ partialtrace(ρ, 2, [4; 4]) ⊗ eye(nE_A) ⊗ eye(nBhat) ⊗ e1 ⊗ partialtrace(ρ, 1, [4; 4]) ⊗ eye(nE_B) * C_AbarBbar) ≥ 0]

  # constraints from purity
  problem.constraints += [trace(C_Abar * W_Abar) ≤ t_A]
  problem.constraints += [trace(C_Abar * W_Bbar) ≤ t_B]

  solve!(problem, SCSSolver(verbose = verbose))

  problem.optval *= nA * nB

  p_succ = trace(eye(nAhat) ⊗ e1 ⊗ partialtrace(ρ, 2, [4; 4]) ⊗ eye(nE_A) ⊗ eye(nBhat) ⊗ e1 ⊗ partialtrace(ρ, 1, [4; 4]) ⊗ eye(nE_B) * C_AbarBbar.value)
  F = problem.optval / p_succ

  return (C_AbarBbar.value, problem, F, p_succ)
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
  problem.constraints += ([0 ⪯ W_Xbar])
  problem.constraints += ([W_Xbar ⪯ eye(size(W_Xbar)...)])
  problem.constraints += ([trace(W_Xbar) == size(C_X)[1] - 1])

  solve!(problem, SCSSolver(verbose = verbose))
  return (W_Xbar.value, problem.optval)
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
