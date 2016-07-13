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

function RainsProbWithPT(ρ, nA, nB, K, δ_min, δ_max, warmstart = false)

	# Check whether ρ is a quantum state
	if isQuantumState(ρ) == 0
		error("The input state must be a valid quantum state.")
	end

	# Check whether dimensions match
	(d,db) = size(ρ)
	if (d != nA*nB)
		error("Input state doesn't match given dimensions." , d , "≠", nA*nB)
	end

	# define the identity matrix on the whole space
	id = eye(d,d)

	ρ_transposed = ρ'
	ρ = ρ

	# define the variable F and the one which is going to be the
	# partial transpose
	# D = Variable(2d, 2d)
	# DPT = Variable(2d, 2d)
	# E = Variable(2d, 2d)
	# EPT = Variable(2d, 2d)
	D = Semidefinite(d)
	E = Semidefinite(d)

	DPT = Semidefinite(d)

	# define the objective
	println("Defining objective...")
	problem = maximize(nA * nB * trace(D * ρ_transposed))
	println("Defined objective.")
	println()

	problem.constraints += ([(id/d - (D+E)) in :SDP])

	# Constraints relating to the PPT Condition
	println("Defining partial transpose constraints...")

	PTM = PTmatrix((d, d), (nB, nB))

	problem.constraints += ([DPT + reshape(PTM*vec(E), d, d)/(K+1)) in :SDP])
	println("...")
	problem.constraints += ([(- reshape(PTM*vec(D), d, d) + reshape(PTM*vec(E), d, d)/(K-1)) in :SDP])
	println("Defined partial transpose constraints.")

	# Constraint coming from the success probability
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E)) ≤ δ_max])
	problem.constraints += ([nA * nB * trace(ρ_transposed*(D+E)) ≥ δ_min])

	println("")
	println("Started solver")

	solve!(problem, SCSSolver(verbose = false))

	println("Finished solver")
	p_succ = nA * nB * trace(ρ_transposed*(D.value+E.value))
	F = problem.optval / p_succ
	return (F, p_succ)
end

# p = 0.5
# nCopies = 2
# rho = sortAB(copies(ronaldState(p),  nCopies), 2, nCopies)
# X = RainsProbRailWithPT(rho, 2^nCopies, 2^nCopies, 2, 0.45, 0.55)
# print(X)
