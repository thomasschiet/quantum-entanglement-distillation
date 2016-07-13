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

function RainsProb(ρ, nA, nB, K, δ_min, δ_max, warmstart = false)

	# Check whether ρ is a quantum state
	if isQuantumState(ρ) == 0
		error("The input state must be a valid quantum state.");
		return;
	end

	# Check whether dimensions match
	(d,db) = size(ρ);
	if (d != nA*nB)
		error("Input state doesn't match given dimensions." , d , "≠", nA*nB);
		return;
	end

	# define the identity matrix on the whole space
	id = complexToReal(eye(d,d));

	ρ_transposed = complexToReal(ρ')
	ρ = complexToReal(ρ);

	# define the variable F and the one which is going to be the
	# partial transpose
	# D = Variable(2d, 2d);
	# DPT = Variable(2d, 2d);
	# E = Variable(2d, 2d);
	# EPT = Variable(2d, 2d);
	D = Variable(2d, 2d);
	DPT = Variable(2d, 2d);
	E = Variable(2d, 2d);
	EPT = Variable(2d, 2d);

	# define the objective
	problem = maximize(nA * nB * trace(D * ρ_transposed) / 2);

	problem.constraints += ([D in :SDP]);
	problem.constraints += ([E in :SDP]);
	problem.constraints += ([(id/d - (D+E)) in :SDP]);

	# give constrains to the used isomorphism and to the matrix being hermitian
	for r = 1:d
		for c = 1:d
			# Define coordinates for the relative blocks
			# r is row, c is column
			# NE | NW     Re | -Im
			# -------  =  --------
			# SE | SW     Im | Re
			rNE = r;
			cNE = c

			rNW = r;
			cNW = c + d

			rSE = r + d;
			cSE = c

			rSW = r + d;
			cSW = c + d

			# NE and SW should be the same because they are both the real part of F
			problem.constraints += ([D[rNE, cNE] ==  D[rSW, cSW]])
			# NW and SE should be
			problem.constraints += ([D[rNW, cNW] == -D[rSE, cSE]])

			# do the same for matrix E
			problem.constraints += ([E[rNE, cNE] ==  E[rSW, cSW]])
			problem.constraints += ([E[rNW, cNW] == -E[rNW, cNW]])

			# The matrix should also be Hermitian.
			# Inequality is because of symmetry around
			# the diagonal of the matrix.
			if r ≥ c
				problem.constraints += ([D[rNE, cNE] ==  D[cNE, rNE]])
				problem.constraints += ([D[rNW, cNW] == -D[cNW, rNW]])
				problem.constraints += ([E[rNE, cNE] ==  E[cNE, rNE]])
				problem.constraints += ([E[rNW, cNW] == -E[cNW, rNW]])
			end
		end
	end

	# Define DPT to be the partial transpose of D
	# and EPT to be the partial transpose of E
	for r = 0:(nA-1)
	for s = 0:(nA-1)
		offsetRow = r * nA;
		offsetCol = s * nA;
		for k = 1:nB
		for l = 1:nB
			p1 = offsetRow + k;
			p2 = offsetCol + l;
			q1 = offsetRow + l;
			q2 = offsetCol + k;
			problem.constraints += ([DPT[q1,q2] == D[p1,p2]]);
			problem.constraints += ([DPT[p1,p2] == D[q1,q2]]);
			problem.constraints += ([EPT[q1,q2] == E[p1,p2]]);
			problem.constraints += ([EPT[p1,p2] == E[q1,q2]]);
			problem.constraints += ([DPT[q1+d,q2+d] == D[p1+d,p2+d]]);
			problem.constraints += ([DPT[p1+d,p2+d] == D[q1+d,q2+d]]);
			problem.constraints += ([EPT[q1+d,q2+d] == E[p1+d,p2+d]]);
			problem.constraints += ([EPT[p1+d,p2+d] == E[q1+d,q2+d]]);

			# The same must hold for the imaginary part
			problem.constraints += ([DPT[q1+d,q2] == D[p1+d,p2]]);
			problem.constraints += ([DPT[p1+d,p2] == D[q1+d,q2]]);
			problem.constraints += ([EPT[q1+d,q2] == E[p1+d,p2]]);
			problem.constraints += ([EPT[p1+d,p2] == E[q1+d,q2]]);
			problem.constraints += ([DPT[q1,q2+d] == D[p1,p2+d]]);
			problem.constraints += ([DPT[p1,p2+d] == D[q1,q2+d]]);
			problem.constraints += ([EPT[q1,q2+d] == E[p1,p2+d]]);
			problem.constraints += ([EPT[p1,p2+d] == E[q1,q2+d]]);
		end
		end
	end
	end

	# Constraints relating to the PPT Condition
	problem.constraints += ([(DPT + EPT/(K+1)) in :SDP]);
	problem.constraints += ([(- DPT + EPT/(K-1)) in :SDP]);

	# Constraint coming from the success probability
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E))/2 <= δ_max]);
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E))/2 >= δ_min]);


	solve!(problem, SCSSolver(verbose=false));

	p_succ = nA * nB * trace(ρ_transposed*(D.value+E.value))/2
	F = problem.optval / p_succ;

	F_out = F
	F_in = eprFidelity(ρ)
	if F_out < F_in
		warn("The output fidelity ", F, " is lower than the input fidelity ", )
	end

	if F_in < 0.5 && F_out > 0.5
		warn("The output fidelity ", F_out, " is higher than 0.5 while the input fidelity ", F_in, " is lower than 0.5")
	end

	return (F, p_succ)
end
