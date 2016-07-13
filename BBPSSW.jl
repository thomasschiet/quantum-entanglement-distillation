include("BEPQuantum.jl")
using BEPQuantum

F_input = collect(0.5:0.02:1)
F_output(F) =  (F.^2+((1-F)/3).^2)./(F.^2+2F.*(1-F)/3 + 5((1-F)/3).^2)
p_succ(F) =  (F^2+2F*(1-F)/3 + 5((1-F)/3)^2)

# F(p_succ) = 1/4 .*(1 .+3 .* sqrt(2 .* F .- 1))

x = collect(0.6:0.01:1)
y = 1/4 .* (1 .+ 3 .* sqrt(2 .* x .-1))
y = F_output(y)

i = 1
for asdf in y
  println("(", x[i], ", ", y[i], ")")
  i = i + 1
end

plot(x=x, y=F_output(y), Geom.line)
1/4 .* (1 .+ 3 .* sqrt(collect(.6:.1:1)))
F_res = F(collect(.6:.1:1))

plot(x=0.6:0.1:1, y=F(collect(0.6:0.1:1)))
# F(p) = 1/4 .+ 3/4 .*p
# p(F) = 4/3 .*(F - 1/4)


F_output(F_input)
ps = p(F_input)

Fs = []

for p_input in ps
  δ_min = 0
  δ_max = p_succ(F(p_input))
  ρ = wernerState(p_input)
  state = sortAB(copies(ρ, 2), 2, 2);

  print("computing (δ_max, δ_min, p, state) = (", δ_max, ", ", δ_min, ", ", p, ", Double Werner)")

  (F_out, p_succ_out) = RainsProb(state, 4, 4, 2, δ_min, δ_max)
  Fs = [Fs, F_out]

  println("found (F_out, p_succ_out) = (", F_out, ", ", p_succ_out, ")")
end

using Gadfly

plot(x=F_input, y = Fs, Geom.line)

i = 1
for x in Fs
  println("(", F_input[i], ", ", x, ")")
  i = i + 1
end
