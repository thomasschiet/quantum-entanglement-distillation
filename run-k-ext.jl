include("BEPQuantum.jl")

using BEPQuantum

n = 2
state = wernerState(0.6)
PPTprogrammeNoTwirling1Ext(copies(state, n), 2n, 2n, 2, 2, 0.6, verbose = true)
