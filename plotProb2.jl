include("quantum.jl")

using Gadfly
using DataFrames

# Range over wernerstates with p between 0 and 1
ps = 0:0.1:1
F = zeros(ps)

K = 2
δs = 0.1:0.1:0.9

# Make the actual states
ρ = map(wernerState, ps)
# Do the optimisation
Fs(δ) = map((state) -> RainsProb(state, 2, 2, K, δ)[3].optval, ρ)
df = DataFrame(p = ps)

df[:F2] = Fs(0.01)
df[:F3] = Fs(0.2)
df[:F4] = Fs(0.26)
df[:F5] = Fs(0.99)

df

# Plot Fidelity vs p and save it
FidelityvsP = plot(df,
  Coord.cartesian(xmin=0, xmax=1, ymin=0, ymax=1),
  layer(x="p", y="F2", Geom.line),
  layer(x="p", y="F3", Geom.line),
  layer(x="p", y="F4", Geom.line),
  layer(x="p", y="F5", Geom.line),
  Theme (
    panel_fill=color("#FFFFFF")
  ),
  Stat.xticks(ticks=[0:0.1:1]),
  Stat.yticks(ticks=[0:0.1:1]),
  Guide.XLabel("p"),
  Guide.YLabel("Fidelity"),
  Guide.Title("Fidelity of Werner states (K ranges from 1 to 5)")
)
draw(PNG("plots/wernerstates.png", 4inch, 3inch), FidelityvsP)
