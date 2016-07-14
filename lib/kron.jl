using Convex

import Base.kron

export kron

function kron(a::Convex.Constant, b::Convex.Variable)
  Rs = AbstractExpr[]
  for i in 1:size(a)[1]
    Vs = Convex.AbstractExpr[]
    for j in 1:size(a)[2]
      push!(Vs, a[i, j] * b)
    end
    push!(Rs, foldl(hcat, Vs))
  end
  return foldl(vcat, Rs)
end
