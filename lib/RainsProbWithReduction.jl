#
# Implements the Rains bound for distillable entanglement but now
# with an allowed failure probability (which in turn may give higher fidelity)
#
# Also implements reduction criterion
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

export RainsProbWithReduction

function RainsProbWithReduction(ρ, nA, nB, K, δ_min, δ_max; verbose = false)

	# Check whether ρ is a quantum state
	@assert isQuantumState(ρ) "The input state must be a valid quantum state."

	# Check whether dimensions match
	(d, db) = size(ρ)
	@assert d == nA * nB "Input state doesn't match given dimensions." , d , "≠", nA*nB

	# define the identity matrix on the whole space
	id = eye(d)

	# define the variable F and the one which is going to be the
	# partial transpose
	D = Semidefinite(d)
	E = Semidefinite(d)

	EPT = ptranspose(E)
	DPT = ptranspose(D)

	ρEpr = epr * epr';
	C = ρEpr ⊗ E + (eye(nA) - ρEpr)/(K^2 - 1) ⊗ E

	# define the objective
	problem = maximize(trace(nA * nB * D * ρ'));

	problem.constraints += (id/d - (D+E)) in :SDP

	# Constraints relating to the PPT Condition
	problem.constraints += (DPT + EPT/(K+1)) in :SDP
	problem.constraints += (- DPT + EPT/(K-1)) in :SDP

	# Constraint coming from the success probability
	problem.constraints += nA * nB * trace(ρ'*(D+E)) ≤ δ_max
	problem.constraints += nA * nB * trace(ρ'*(D+E)) ≥ δ_min

	problem.constraints +=             partialtrace(C, 1, [nA * 2; nB * 2]) ⊗ eye(nB*2) ⪰ C
	problem.constraints += eye(nA*2) ⊗ partialtrace(C, 2, [nA * 2; nB * 2])             ⪰ C

	solve!(problem, SCSSolver(verbose = verbose, eps = 5e-3, max_iters = 1e7))

	p_succ = nA * nB * trace(ρ'*(D.value + E.value))
	F = problem.optval / p_succ

	return (problem, F, p_succ, D.value, E.value)
end
