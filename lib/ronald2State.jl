#
# Outputs a ronald state of the form p EPR + (1-p) max mixed
#
# Inputs: p

function ronald2State(p, φ = 0)

	e0 = eVec(2, 1)
	e1 = eVec(2, 2)
	# The EPR Pair
	epr = (e0 ⊗ e1 + e^(im*φ) * e1 ⊗ e0 )/√2;

	v00 = e0 ⊗ e1;

	out = p * epr*epr' + (1-p) * v00*v00';

  return out

end
