#
# Checks whether the input is a valid unitary
#

function isUnitary(U)

	# defined the cutoff when we consider something to be positive
	epsilon = 10.0^(-15);

	# Check dimensions
	(d,da) = size(U);
	if (d != da) 
		error("Input is not a square matrix");
		return false;
	end

	# Check unitarity
	if (U*U' != eye(d))
		return false
	end
	if (U'*U != eye(d))
		return false
	end

	return true;
end

	
