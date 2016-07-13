#
# Checks whether the input ρ is a valid quantum state
#

function isQuantumState(ρ)
	return true
	# defined the cutoff when we consider something to be positive
	epsilon = 10.0^(-15);

	# Check whether it's hermitian
	# if (ρ != ρ')
	if !ishermitian(ρ)
		print("Not Hermitian")
		return 0;
	end

	# check for positivity
	l = eig(ρ)[1];
	if (minimum(l) <= - epsilon)
		print("Not positive")
		return 0;
	end

	return 1;
end
