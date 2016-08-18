#
# Produces n copies of a density matrix
#
# Input:
# ρ	density matrix
# n	number of copies we want
#
export copies

function copies(ρ, n)
	return reduce(⊗, repeated(ρ, n))
end
