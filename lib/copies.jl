#
# Produces n copies of a density matrix
#
# Input:
# ρ	density matrix
# n	number of copies we want
#
export copies

function copies(ρ::AbstractMatrix, n::Int)
	return reduce(⊗, repeated(ρ, n))
end

function copies(ρ::AbstractMatrix, n::Float64)
	return copies(ρ, convert(Int, n))
end
