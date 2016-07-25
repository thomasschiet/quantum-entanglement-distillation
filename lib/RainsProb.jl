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
# (F, p_succ) where:
# F	fidelity
# p_succ change of success

using Convex
using SCS

export RainsProb

function RainsProb(ρ, nA, nB, K, δ_min, δ_max; verbose = false)

	# Check whether ρ is a quantum state
	if isQuantumState(ρ) == 0
		error("The input state must be a valid quantum state.");
		return;
	end

	# Check whether dimensions match
	(d, db) = size(ρ);
	if (d != nA*nB)
		error("Input state doesn't match given dimensions." , d , "≠", nA*nB);
		return;
	end

	# define the identity matrix on the whole space
	id = eye(d)

	# define the variable F and the one which is going to be the
	# partial transpose
	M = Semidefinite(d)
	MPT = Semidefinite(d)

	E = Semidefinite(d)
	EPT = Semidefinite(d)

	p_succ = nA * nB * trace(ρ'*(M+E))

	# define the objective
	problem = maximize(nA * nB * trace(ρ' * M))


	EPT = partialtranspose(E)
	MPT = partialtranspose(M)

	problem.constraints += ([M + E ≤ id/d])

	# Constraints relating to the PPT Condition
	problem.constraints += ([MPT + EPT/(K+1) ≥ 0])
	problem.constraints += ([- MPT + EPT/(K-1) ≥ 0])

	# Constraint coming from the success probability
	problem.constraints += ([p_succ <=  δ_max])
	problem.constraints += ([p_succ >=  δ_min])

	solve!(problem, SCSSolver(verbose = verbose));

	p_succ = nA * nB * trace(ρ'*(M.value + E.value))

	F = problem.optval / p_succ;

	return (problem, F, p_succ)
end
