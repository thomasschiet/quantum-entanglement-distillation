include("BEPQuantum.jl")

using BEPQuantum

ρ = sortAB(copies(wernerState(0.7), 2), 2, 2)
purity(ρ)
