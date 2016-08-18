include("BEPQuantum.jl")

using BEPQuantum

n = 2
state = ronaldState(0.4)
δ = 0.5799999999999995

println(deutsch(state))
bennett(full(state))

# 0.8008173307612931
# 0.580046288135209
(problem, F, p_succ, Wa) = PPTprogrammeNoTwirling1ExtPermSym(copies(state, n), 2n, 2n, 2, 2, δ, verbose = true)
println(F)
println(p_succ)
println(Wa.value)
println(maximum(Wa.value))
println(sum(abs(Wa.value)))

# 0.800798757626111
# 0.5799992958755701
(problem, F, p_succ) = PPTprogrammeNoTwirling1ExtPermSymOnlySym(copies(state, n), 2n, 2n, 2, 2, δ, verbose = true)
println(F)
println(p_succ)

sortAB(copies(state, n), 2, n)
(problem, F, p_succ) = RainsProb(sortAB(copies(state, n), 2, n), 2n, 2n, 2, 0, δ, verbose = true)
println(F)
println(p_succ)
