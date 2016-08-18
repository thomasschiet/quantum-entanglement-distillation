export random_unitary, random_densitymatrix

function random_unitary(dim::Integer; real::Bool = true, complex::Bool = false)
  @assert dim ≥ 1

  generated = randn(dim, dim)

  if complex || !real
    generated += im * randn(dim, dim)
  end

  (Q, R) = qr(generated)

  F = diagm(diag(R))

  R = sign(diagm(diag(R)))
  R[find((x)->  x == 0, R)] = 1

  F = diagm(diag(R))

  return Q*F
end

function random_densitymatrix(dim::Integer; real::Bool = true, complex::Bool = false)
  @assert dim ≥ 1

  U = random_unitary(dim, real = real, complex = complex)
  U *= U'
  return U/trace(U)

end
