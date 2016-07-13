module BEPQuantum
	export ⊗, ⊕, issquare
	⊗(x, y) = kron(x, y)
	⊕(x::Matrix, y::Matrix) = [x zeros(size(x)[1], size(y)[2]); zeros(size(y)[1], size(x)[2])]
	issquare(x::AbstractArray) = size(x)[1] == size(x)[2]

	include("lib/eVec.jl");
	include("lib/commonStates.jl");
	include("lib/bellState.jl");
	include("lib/braPT.jl");
	include("lib/ketPT.jl");
	include("lib/isQuantumState.jl");
	include("lib/isPPT.jl");
	include("lib/eprFidelity.jl");

	include("lib/Rains.jl");
	include("lib/RainsProb.jl");

	include("lib/copies.jl");
	include("lib/isUnitary.jl");
	include("lib/sortAB.jl");

	include("lib/realToComplex.jl");
	include("lib/complexToReal.jl");

	include("lib/partialtrace.jl");
	include("lib/partialtranspose.jl");
	include("lib/permutesystems.jl");

	include("lib/constructEpsnet.jl");
	include("lib/purity.jl");

	include("lib/colors.jl");
end
