#
# Produces n copies of a density matrix
#
# Input:
# ρ	density matrix
# n	number of copies we want
#
export copies

function copies(ρ,n)

	out = copy(ρ);
	for k=1:(n-1)
		out = kron(out,ρ);
	end

	return out;
end
