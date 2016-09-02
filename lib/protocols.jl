export bennett
export deutsch
export epl

import Base.convert
convert(::Type{Float64}, x::Array{Float64,1}) = x[1]
convert(::Type{Float64}, x::Array{Float64,2}) = x[1]

toMatrix(x::VecOrMat) = x*x'

e0 = [1; 0]
e1 = [0; 1]

⊗(x, y) = kron(x, y)
Φ_plus = 1/√2 * (e0 ⊗ e0 + e1 ⊗ e1)
Φ_min  = 1/√2 * (e0 ⊗ e0 - e1 ⊗ e1)
Ψ_plus = 1/√2 * (e0 ⊗ e1 + e1 ⊗ e0)
Ψ_min  = 1/√2 * (e0 ⊗ e1 - e1 ⊗ e0)
basis = Any[Φ_plus, Φ_min, Ψ_plus, Ψ_min]

# twirl
function toBellBasis(ρ::AbstractMatrix)
  output = zeros(size(ρ))
  i = 1
  for x in basis
    j = 1
    for y in basis
      output[i, j] = (x' * ρ * y)[1]
      j += 1
    end
    i += 1
  end
  return output
end

function deutsch(ρ::AbstractMatrix)
  if fidelity(rot(ρ)) > fidelity(ρ)
    ρ = rot(ρ)
  end
  M = toBellBasis(ρ)
  p = (M[1, 1] + M[4, 4])^2 + (M[2, 2] + M[3, 3])^2
  F = (M[1, 1]^2 + M[4, 4]^2)/p
  return (F, p)
end


∫(x::AbstractMatrix, y::Vector) = twirl(x, y)
function twirl(ρ::AbstractMatrix, twirl::Vector)
  φ = zeros(size(ρ))
  for op in twirl
    φ += op * ρ * op'
  end
  φ /= length(twirl)
  return φ
end

function polarize(ρ::AbstractMatrix)
  # X Pauli
  X = [0 1; 1 0]
  # Z Pauli
  Z = [1 0; 0 -1]

  id = eye(4)
  K1 = X ⊗ X
  K2 = Z ⊗ Z
  K1K2 = K1*K2

  # Perform the twirl
  twirlops = Any[id, K1, K2, K1K2]
  return ∫(ρ, twirlops)
end

# Calculate fidelity of a matrix ρ
fidelity(ρ::AbstractMatrix) = Float64(real(basis[1]' * ρ * basis[1]))

function bennett(F::Real)
  p_succ = F^2 + 2F*(1-F)/3 + 5((1-F)/3)^2
  F_out = (F^2 + ((1-F)/3)^2)/p_succ
  return (F_out, p_succ)
end

function bennett(ρ::AbstractMatrix)
  if fidelity(rot(ρ)) > fidelity(ρ)
    ρ = rot(ρ)
  end
  ρ = polarize(ρ)

  F = fidelity(ρ)

  return bennett(F)
end

function rot(ρ::AbstractMatrix)
  # X Pauli
  X = [0 1; 1 0]
  # Z Pauli
  Z = [1 0; 0 -1]
  id = eye(2)
  return (id ⊗ X) * ρ * (id ⊗ X)'
end

function epl(p::Number)
  return (1, 0.5*p^2)
end
