#
# Given ρ_AB^tensor n, we want to shuffle the A and B systems
# such that we get ρ_A1...An,B1...Bn, ie all the A's and B's have been
# sorted to both sides.
#
# Inputs:
# ρBig	the large state, doesnt have to be of tensor form
#		but A's and B's are assumed to be equally large everywhere
# d		dimension of A, which is the same as the dimension of B
# n		number of "copies"
#
#
# Outputs:
# ρSorted	a new state in which all A's and B's are sorted towards
#		the sides
#
export sortAB

function sortAB(ρBig, d, n)

	# Check whether the input is a valid state and the dimensions match
	if (isQuantumState(ρBig) == 0)
		error("The input is not a valid quantum state.");
		return;
	end
	(d1,d2) = size(ρBig);
	if (d1 != (d^2)^n)
		error("Input dimensions don't match." , d1, "≠", (d^2)^n);
		return;
	end

	half = floor(n/2);

	# Build the unitary that shuffles the systems around
	# depending on wether n is even or odd, we need an extra
	# identity term in the middle

	U = zeros((d^2)^n,(d^2)^n);
	if iseven(n)
		for j=1:d
		for k=1:d
			# Swap j with k
			jk = eVec(d,j) * eVec(d,k)';
			kj = eVec(d,k) * eVec(d,j)';
			base1 = kron(eye(d),jk);
			base2 = kron(kj,eye(d));

			U = U + kron(copies(base1,half),copies(base2,half));
		end
		end
	else
		for j=1:d
		for k=1:d
			# Swap j with k
			jk = eVec(d,j) * eVec(d,k)';
			kj = eVec(d,k) * eVec(d,j)';
			base1 = kron(eye(d),jk);
			base2 = kron(kj,eye(d));

			U = U + kron(copies(base1,half),kron(eye(d^2),copies(base2,half)));
		end
		end
	end


	ρOut = copy(ρBig);
	ρOut = U*ρOut*U';

	return ρOut;
end
