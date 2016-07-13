#
# Outputs a ronald state of the form p EPR + (1-p) max mixed
#
# Inputs: p

function ronaldState(p)

	e0 = eVec(2, 1)
	e1 = eVec(2, 2)

	# The EPR Pair
	epr = (e0 ⊗ e0 + e1 ⊗ e1)/√2;
	v11 = e0 ⊗ e0;
	out = p * epr*epr' + (1-p) * v11*v11';

	println(out)

	return out;

end
