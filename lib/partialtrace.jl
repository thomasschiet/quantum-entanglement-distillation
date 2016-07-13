#
# Partial trace
#
# Inputs:
# ρ			Square matrix
# dims	array with dimensions of systems
# sys		system to trace out
#
# Outputs:
# The partial trace of ρ where system `sys` was traced out
#

# include("copies.jl")
# ⊗(x, y) = kron(x, y)
using Convex
using SCS

import Base.sign
import Convex.monotonicity
import Convex.conic_form!
import Convex.curvature

export ptrace, partialtrace

type PartialTraceAtom <: AbstractExpr
	head::Symbol
	id_hash::UInt64
	children::Tuple{AbstractExpr}
	size::Tuple{Integer, Integer}
	sys::Integer
	dims::Vector

	function PartialTraceAtom(x::AbstractExpr, sys::Integer, dims::Vector)
		if x.size[1] ≠ x.size[2]
			error("Only square matrices are supported")
		end

		if ! (1 ≤ sys ≤ length(dims))
			error("Invalid system, should between 1 and ", length(dims))
		end

		if x.size[1] ≠ prod(dims)
			error("Dimension of system doesnt correspond to dimension of subsystems")
		end

		children = (x, )
		newsize = (round(Int, x.size[2]/dims[sys]), round(Int, x.size[1]/dims[sys]))

		return new(:partialTrace,
		hash(children),
		children,
		newsize,
		sys,
		dims)
	end
end

function sign(x::PartialTraceAtom)
	 return NoSign()
	 return sign(x.children[0])
 end

function curvature(x::PartialTraceAtom)
  return ConstVexity()
end

function monotonicity(x::PartialTraceAtom)
	return (NoMonotonicity(),)
end

function conic_form!(x::PartialTraceAtom, unique_conic_forms::UniqueConicForms)
	if !has_conic_form(unique_conic_forms, x)
		sys = x.sys
		dims = x.dims
		function entry(ρ, j::Integer)
			bra = speye(1)
			ket = speye(1)
			i_sys = 1
			for dim in dims
				if i_sys == sys
					vO = sparsevec([j], [1], dim);
					bra = kron(bra, vO')
					ket = kron(ket, vO)
				else
					bra = kron(bra, speye(dim))
					ket = kron(ket, speye(dim))
				end
				i_sys += 1
			end
			return bra * ρ * ket
		end

		objective = conic_form!(sum([entry(x.children[1], j) for j in 1:dims[sys]]), unique_conic_forms)
		cache_conic_form!(unique_conic_forms, x, objective)
	end
	return get_conic_form(unique_conic_forms, x)
end

# General function, works in Convex.jl and with normal matrices
function partialtrace(x, sys::Integer, dims::Vector)
	function entry(ρ, j::Integer)
		bra = speye(1)
		ket = speye(1)
		i_sys = 1
		for dim in dims
			if i_sys == sys
				vO = sparsevec([j], [1], dim)
				bra = kron(bra, vO')
				ket = kron(ket, vO)
			else
				bra = kron(bra, speye(dim))
				ket = kron(ket, speye(dim))
			end
			i_sys += 1
		end
		return bra * ρ * ket
	end

	return sum([entry(x, j) for j in 1:dims[sys]])
end


# ignore untill it's fixed
# partialtrace(x::AbstractExpr, sys::Integer, dim::Vector) = PartialTraceAtom(x, sys, dim)
# ptrace(x::Variable, sys::Integer, dim::Vector) = PartialTraceAtom(x, sys, dim)

ptrace(x::AbstractExpr, sys::Integer, dim::Vector) = partialtrace(x, sys, dim)
ptrace(x::AbstractArray, sys::Integer, dim::Vector) = partialtrace(x, sys, dim)
