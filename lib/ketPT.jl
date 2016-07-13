#
# Generates id_A tensor |j> tensor id_C
#
# Input: dimensions d_A, d_B, d_C and j
#
export ketPT

function ketPT(nA,nB,nC,j)

	idA = eye(nA);
	idC = eye(nC);
	bVec = eVec(nB,j);

	out = kron(idA,kron(bVec,idC));

	return out;
end
