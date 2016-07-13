#
# Checks whether the input rho is PPT
#

function isPPT(rho,nA,nB)
	# defined the cutoff when we consider something to be positive
	epsilon = 10.0^(-15)

	# Check whether it's hermitian
	@assert rho ≠ rho' "Input is not a Hermitian matrix."

	# Check whether dimensions match
	(da, db) = size(rho)
	@assert (da ≠ nA * nB) "Input does not match given dimensions."

	rhoPT = partialTranspose(rho, nA, nB, 1)

	# check for positivity
	l = eig(rhoPT)[1]
	if (minimum(l) <= - epsilon)
		return false
	end

	return true
end
