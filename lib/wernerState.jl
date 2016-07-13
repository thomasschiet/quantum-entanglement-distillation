#
# Outputs a werner state of the form p EPR + (1-p) max mixed
#
# Inputs: p
# module Quantum
	# export wernerState
	function wernerState(p)

		# The EPR Pair
		epr = (kron(eVec(2,1),eVec(2,1)) + kron(eVec(2,2),eVec(2,2)))/sqrt(2);

		out = p * epr*epr' + (1-p) * eye(4)/4;

		return out;

	end
# end
