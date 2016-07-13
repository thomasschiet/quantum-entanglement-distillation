#
# Make some plots of the fidelities obtained in entanglement distillation
#
using Convex
include("quantum.jl")
include("PTmatrix.jl")

n = 3
state = sortAB(copies(ronaldState(0.6), n), 2, n);
state = state
x = Rains(state, 2^n, 2^n, 2)
print(x.optval)


#
# t = linspace(0,1,11)
#
# δ = 0.1
#
# y = [];
# z = [];
# for j = 1:11
# 	state = sortAB(copies(ronaldState(t[j]),2),2,2);
# 	pr = RainsProb(state, 4, 4, 2, δ)
# 	y = [y; pr];
# 	z = [z; eprFidelity(ronaldState(t[j]))];
# end
#
# using Gadfly
# using DataFrames
#
# RainsProb(sortAB(copies(ronaldState(0.1),2),2,2), 4, 4, 2, δ)
#
# df = DataFrame(p = t)
# df[:y] = y
# df[:z] = z
#
#
# df
#
# FidelityvsP = plot(df,
#   Coord.cartesian(xmin=0, xmax=1, ymin=0, ymax=1),
#   layer(x="p", y="y", Geom.line, Theme(default_color=color("red"))),
#   layer(x="p", y="z", Geom.line, Theme(default_color=color("blue"))),
#   Theme (
#     panel_fill=color("#FFFFFF")
#   ),
#   Stat.xticks(ticks=[0:0.1:1]),
#   Stat.yticks(ticks=[0:0.1:1]),
#   Guide.XLabel("p"),
#   Guide.YLabel("Fidelity"),
#   Guide.Title("Fidelity with and without RainsProb for Rolandstate"))
#
# draw(PNG("plots/ronald_delta_0_9_c_0_5.png", 8inch, 6inch), FidelityvsP)
