include("../BEPQuantum.jl")

# it should work in Convex.jl
using BEPQuantum
using FactCheck
using Convex
using SCS

TOL = 1e-3

include("permutesystems.jl")
include("partialtranspose.jl")
include("partialtrace.jl")
include("kron.jl")
include("directsum.jl")
include("rains.jl")
include("purity.jl")
include("protocols.jl")
