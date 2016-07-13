#
# Computes the fidelity with the EPR pair
#
export eprFidelity

function eprFidelity(ρ)

	# Check ρ is a valid input
	(d,db) = size(ρ)
	if (d != 4) || (isQuantumState(ρ) == 0)
		error("The input must be a 2x2 quantum state");
		return;
	end

	# The EPR Pair
	epr = (kron(eVec(2,1),eVec(2,1)) + kron(eVec(2,2),eVec(2,2)))/sqrt(2);

	out = (epr' * ρ * epr)[1];

	return out;

end
