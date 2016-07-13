using Gadfly
using Colors

id = eye(2)
pauli1 = [0 1; 1 0]
pauli2 = [0 -im; im 0]
pauli3 = [1 0; 0 -1]

function f(x::Number, y::Number)
  return max(minimum(eig(kron(id, id) + x*kron(pauli1, pauli1) + y*kron(pauli1, pauli3))[1]), 0);
end

function f(xs::Array, ys::Array)
  res = []
  i = 1
  for x in xs
    y = ys[i]
    res = [res; f(x, y)]
    i += 1
  end
  return res
end

f_max(t) = maximum(eig(id + t*pauli3 + (t-1)*pauli2)[1])
f_min(t) = minimum(eig(id + t*pauli3 + (t-1)*pauli2)[1])

x = repeat(collect(-1:0.1:1), inner=[21])
y = repeat(collect(-1:0.1:1), outer=[21])

plot(x=x, y=y, color=f(x, y), Geom.rectbin,
     Scale.color_continuous())
