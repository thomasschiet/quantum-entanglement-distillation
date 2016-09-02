include("BEPQuantum.jl")

using BEPQuantum

# n = 2
# p = 0.4
# state = wernerState(p)
# δ = bennett(state)[2]
# # δ = 0.5
#
# verbose = true
# eps = 1e-4
#
# # RONALDSTATE
# # p = 0.4
# # δ = 0.5
# # 0.8647718648161498
# # 0.8647716307665765
# # 0.8647716650946311
#
# # p = 0.7
# # δ = 0.5
# # 0.976969282039843 Solve time: 4.23e-01s
# # 0.976955153941947 Solve time: 3.35e+01s
# # 0.976959426401371 Solve time: 1.67e+01s
#
# # RONALD2STATE (p = 0.4, δ = 0.5, eps = 1e-7)
# # 0.6377777783763448 Solve time: 4.54e-01s
# # 0.6377777630083299 Solve time: 2.95e+01s
# # 0.6377777814871788 Solve time: 1.99e+01s
#
# # RONALD2STATE (p = 0.7, δ = 0.5, eps = 1e-4)
# # 0.9828446157835343 Solve time: 1.78e+00s
# #
#
# # WERNERSTATE (p = 0.4, δ = 0.5, eps = 1e-7)
# # 0.5990099077656231 Solve time: 2.45e+00s
# # 0.5944352306876611 Solve time: 6.36e+02s
# # 0.594304894578052
# # 0.5944352370371205 Solve time: 2.68e+02s
#
# # WERNERSTATE (p = 0.4, δ = 0.5, eps = 1e-4)
# # 0.5990046426780034 Solve time: 8.67e-01s
# # -
# # 0.5944669555720423 Solve time: 8.64e+01s
#
#
#
# # δ = 0.744964339650292
# # p = 0.7
# # 0.8691833585699156 Solve time: 4.14e+00s
# # 0.869000503378117  Solve time: 2.00e+02s
# # 0.8690089459283302 Solve time: 8.31e+01s
#
#
# 0.8691833585699156 - 0.869000503378117
#
# # RainsProb
# # PPTprogrammeNoTwirling1ExtPermSym
# # PPTprogrammeNoTwirling1ExtPermSymOnlySym
# # RONALDSTATE
# # p = 0.4
# # δ = 0.5
# # 0.8647718648161498
# # 0.5000028114571987
# # eps = 1e-6
# # 0.8647718760173041
# # 0.4999999224956734
# println("RainsProb")
# (problem, F, p_succ, Wa) = RainsProb(state, 2n, 2n, 2, 2, 0, δ, verbose = verbose, eps = eps)
# println(F)
# println(p_succ)
#
# # RONALDSTATE
# # p = 0.4
# # δ = 0.5
# # 0.8646009017141222
# # 0.5001020502434709
# # eps = 1e-6 Solve time: 2.43e+01s
# # 0.8647697483090633
# # 0.5000012566843111
#
# # 0.8008173307612931
# # 0.580046288135209
# # println("PPTprogrammeNoTwirling1ExtPermSym")
# # (problem, F, p_succ, Wa) = PPTprogrammeNoTwirling1ExtPermSym(copies(state, n), 2n, 2n, 2, 2, δ, verbose = verbose, eps = eps)
# # println(F)
# # println(p_succ)
#
# # RONALDSTATE
# # p = 0.4
# # δ = 0.5
# # 0.8645646588280245
# # 0.5001246376929067
#
# # eps = 1e-6
# # 0.8647698317704671
# # 0.5000012165984381
#
# # 0.800798757626111
# # 0.5799992958755701
# println("PPTprogrammeNoTwirling1ExtPermSymOnlySym")
# (problem, F, p_succ) = PPTprogrammeNoTwirling1ExtPermSymOnlySym(copies(state, n), 2n, 2n, 2, 2, δ, verbose = verbose, eps = eps)
# println(F)
# println(p_succ)
#
# # (problem, F, p_succ) = PPTprogrammeNoTwirling2Ext(copies(state, n), 2n, 2n, 2, 2, δ, verbose = verbose, eps = eps)
# # println(F)
# # println(p_succ)
#
# # 0.9222319777245588
# # 0.9134588223911435
#
# 0.9222319777245588-0.9134588223911435

id = eye(2)
X = [0 1; 1 0]
state = ronald2State(0.6)
op = id ⊗ X
eprFidelity(op * state * op)
deutsch(op * state * op)[1]
δ = deutsch(state)[2]
@time x = RainsProb(state, 4, 4, 2, 2, 0, δ)
x[2]
@time x = PPTprogrammeNoTwirling1ExtPermSymOnlySym(copies(state, 2), 4, 4, 2, 2, δ)
x[2]
x
x = PPTprogrammeNoTwirling2Ext(copies(state, 2), 4, 4, 2, 2, δ)
# @time x = RainsProbWithReduction(sortAB(copies(state, 2), 2, 2), 4, 4, 2, 0, δ)
println(x[2])
x[2]
x[1]

# old
x[2]
x[2]
