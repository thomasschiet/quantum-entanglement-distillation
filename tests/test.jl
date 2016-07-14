include("../BEPQuantum.jl")

# it should work in Convex.jl
using BEPQuantum
using FactCheck
using Convex
using SCS

include("permutesystems.jl")
include("partialtranspose.jl")
include("partialtrace.jl")
include("kron.jl")
