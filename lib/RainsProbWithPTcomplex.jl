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

include("quantum.jl")
include("PTmatrix.jl")

function RainsProbRailWithPT(ρ, nA, nB, K, δ_min, δ_max, warmstart = false)

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

	ρ_transposed = complexToReal(ρ');
	ρ = complexToReal(ρ);

	# define the variable F and the one which is going to be the
	# partial transpose
	D = Semidefinite(2d);
	E = Semidefinite(2d);

	# define the objective
	println("Defining objective...")
	problem = maximize(nA * nB * trace(D * ρ_transposed) / 2);
	println("Defined objective.")

	println()
	println("Adding complex constrains")
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
	println("Complex constrains added")

	println()

	problem.constraints += ([(id/d - (D+E)) in :SDP]);

	# Constraints relating to the PPT Condition
	println("Defining partial transpose constraints...")
	problem.constraints += ([(PT(D) + PT(E)/(K+1)) in :SDP]);
	println("...")
	problem.constraints += ([(- PT(D) + PT(E)/(K-1)) in :SDP]);
	println("Defined partial transpose constraints.")

	# Constraint coming from the success probability
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E))/2 ≤ δ_max]);
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E))/2 ≥ δ_min]);

	println("")
	println("Started solver")

	solve!(problem, SCSSolver());

	println("Finished solver")
	p_succ = nA * nB * trace(ρ_transposed*(D.value+E.value))/2
	F = problem.optval / p_succ;

	return (F, p_succ)
end

p = 0.5
nCopies = 3
rho = sortAB(copies(ronaldState(p),  nCopies), 2, nCopies)
X = RainsProbRailWithPT(rho, 2^nCopies, 2^nCopies, 2, 0, 1)
print(X)
