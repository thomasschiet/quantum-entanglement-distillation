#
# Outputs a ronald state of the form p EPR + (1-p) max mixed
#
# Inputs: p
export ronaldStateQutrit

function ronaldStateQutrit(p::Number, φ::Number = 0)

	e0 = eVec(3, 1)
	e1 = eVec(3, 2)
	e2 = eVec(3, 3)
	# The EPR Pair
	epr = (e0 ⊗ e0 + e1 ⊗ e1 + e2 ⊗ e2 )/√3;

	v00 = e0 ⊗ e0;

	out = p * epr*epr' + (1-p) * v00*v00';

  return out

end
