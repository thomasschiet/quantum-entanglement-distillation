#
# Implements the Rains bound for distillable entanglement but now
# with an allowed failure probability (which in turn may give higher fidelity)
#
# Inputs:
# ρ quantum state
# nA  dimension Alice
# nB  dimension Bob (assumed to be the same for now)
# K   output dimension
# δ_min minimum allowed failure probability
# δ_max maximum allowed failure probability
#
# Outputs:
# (problem, F, p_succ)

using Convex
using SCS

export RainsProb

function RainsProb(ρ::AbstractMatrix, nA::Int, nB::Int, n::Int, K::Int, δ_min::Number, δ_max::Number; verbose::Bool = false, eps::Number = 1e-4)

	# copy ρ n times
	ρ = copies(ρ, n)

	# Check whether ρ is a quantum state
	@assert isQuantumState(ρ) "The input state must be a valid quantum state."

	# Check whether dimensions match
	(d, db) = size(ρ)
	@assert d == nA * nB "Input state doesn't match given dimensions." , d , "≠", nA*nB

	# define the identity matrix on the whole space
	id = eye(d)

  #sort all the entries, such that all the systems on A are first and all
  #the systems on B are second
  if n > 1
    rhoSorted = sortAB(ρ, 2, n)
  else
    rhoSorted = ρ
  end

	# define the variable F and the one which is going to be the
	# partial transpose
	D = Semidefinite(d)
	E = Semidefinite(d)

	EPT = ptranspose(E)
	DPT = ptranspose(D)

	# define the objective
	problem = maximize(trace(nA * nB * D * rhoSorted'));
	problem.constraints += (id/d - (D+E)) in :SDP

	# Constraints relating to the PPT Condition
	problem.constraints += (DPT + EPT/(K+1)) in :SDP
	problem.constraints += (- DPT + EPT/(K-1)) in :SDP

	# Constraint coming from the success probability
	problem.constraints += nA * nB * trace(rhoSorted'*(D+E)) ≤ δ_max
	problem.constraints += nA * nB * trace(rhoSorted'*(D+E)) ≥ δ_min

	solve!(problem, SCSSolver(verbose = verbose, eps = eps))
	p_succ = nA * nB * trace(rhoSorted'*(D.value + E.value))
	F = problem.optval / p_succ

	return (problem, F, p_succ, D.value, E.value)
end
