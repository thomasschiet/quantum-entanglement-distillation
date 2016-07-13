# Makes a vector in the standard basis
#
export eVec

function eVec(dim::Int, j::Int)
	@assert dim ≥ 0
	@assert j ≥ 0
	@assert dim ≥ j
	return sparsevec([j], [1], dim)
end
