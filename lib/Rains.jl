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

function Rains(rho, nA::Int, nB::Int, K::Int)

	# Check whether rho is a quantum state
	@assert isQuantumState(rho) "Rho is not a quantum state"

	# Check whether dimensions match
	(d, db) = size(rho);
	@assert d == nA*nB "Input state doesn't match given dimensions.", d, "≠", nA*nB

	# define the identity matrix on the whole space
	id = eye(d,d);

	# define the variable F and the one which is going to be the
	# partial transpose
	F = Semidefinite(d);
	# FPT = Variable(d,d);

	# define the objective
	problem = maximize(trace(F * rho));

	# problem.constraints += ([F in :SDP]);
	problem.constraints += ([(id - F) in :SDP]);

	# Shuffle the entries around
	# for r = 0:(nA-1)
	# for s = 0:(nA-1)
	# 	offsetRow = r * nA;
	# 	offsetCol = s * nA;
	# 	for k = 1:nB
	# 	for l = 1:nB
	# 		p1 = offsetRow + k;
	# 		p2 = offsetCol + l;
	# 		q1 = offsetRow + l;
	# 		q2 = offsetCol + k;
	# 		problem.constraints += ([FPT[q1,q2] == F[p1,p2]]);
	# 		problem.constraints += ([FPT[p1,p2] == F[q1,q2]]);
	# 	end
	# 	end
	# end
	# end

	problem.constraints += ([(PT(F, (nB, nB)) + id/K) in :SDP]);
	problem.constraints += ([(id/K - PT(F, (nB, nB))) in :SDP]);

	solve!(problem,SCSSolver(verbose=true))

	problem.optval

	return problem
end
