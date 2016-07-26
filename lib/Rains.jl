#
# Implements the Rains bound for distillable entanglement
#
# Inputs:
# rho quantum state
# nA  dimension Alice
# nB  dimension Bob (assumed to be the same for now)
# K   output dimension
#
# Outputs:
# problem problem object given by Convex. Optimal value can be obtained
#         by problem.optval
export Rains

using Convex
using SCS

function Rains(rho, nA::Int, nB::Int, K::Int; verbose = true)

	# Check whether rho is a quantum state
	@assert isQuantumState(rho) "Rho is not a quantum state"

	# Check whether dimensions match
	(d, db) = size(rho);
	@assert d == nA*nB "Input state doesn't match given dimensions.", d, "≠", nA*nB

	# define the identity matrix on the whole space
	id = eye(d,d);

	# define the variable F
	F = Semidefinite(d);

	# define the objective
	problem = maximize(trace(F * rho));

	problem.constraints += F ⪯ id;
	problem.constraints += ptranspose(F) + id/K ⪰ 0;
	problem.constraints += id/K - ptranspose(F) ⪰ 0;

	solve!(problem,SCSSolver(verbose = verbose))

	return (problem, F)
end
