include("BEPQuantum.jl")

using BEPQuantum

δ_min = 0
δ_max = 0.5

parameter = 0.5
n_copies = 2

ρ = full(ronald2State(parameter))

reduce(⊗, repeated(ρ, 2))

## Ronald state p = 0.5, δ ∈ [0, 0.5]
# 0.7629379836769585
# 0.5000000273275189
state = sortAB(copies(full(ronald2State(parameter)), n_copies),  2, n_copies)
@time (problem, F, p_succ, D, E) = RainsProb(state, 2^n_copies, 2^n_copies, 2, δ_min, δ_max, verbose = true)
println(F)
println(p_succ)


## Ronald state p = 0.5, δ ∈ [0, 0.5]
# eps = 5e-3
# above one
# 0.7628895488245799
# 0.49999999886237434
#
# lower one
#
#
#
# both
# 0.7628690762884807
# 0.5000000422140645
@time (problem, F, p_succ, D, E) = RainsProbWithReduction(state, 2^n_copies, 2^n_copies, 2, δ_min, δ_max, verbose = true)
println(F)
println(p_succ)
