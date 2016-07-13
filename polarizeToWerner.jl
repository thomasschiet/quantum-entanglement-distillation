# The polarization described here
# is
include("wernerState.jl")
include("ronaldState.jl")
include("ronald2State.jl")
include("eVec.jl")

using Quantum
using Gadfly

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
function toBellBasis(ρ::Matrix)
  bell = zeros(size(ρ))
  i = 1
  for x in basis
    j = 1
    for y in basis
      bell[i, j] = x * ρ * y'
      j += 1
    end
    i += 1
  end
  return bell
end

∫(x::Matrix, y::Vector) = twirl(x, y)
function twirl(ρ::Matrix, twirl::Vector)
  φ = zeros(size(ρ))
  for op in twirl
    φ += op * ρ * op'
  end
  φ /= length(twirl)
  return φ
end

function polarize(ρ::Matrix)
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
fidelity(ρ::Matrix) = Float64(real(basis[1]' * ρ * basis[1]))

function bennett(F::Real)
  p_succ = F^2 + 2F*(1-F)/3 + 5((1-F)/3)^2
  F_out = (F^2 + ((1-F)/3)^2)/p_succ
  return (F_out, p_succ)
end

function bennett(ρ::Matrix)
  if fidelity(rot(ρ)) > fidelity(ρ)
    ρ = rot(ρ)
  end
  ρ = polarize(ρ)

  F = fidelity(ρ)

  return bennett(F)
end

function rot(rho::Matrix)
  # X Pauli
  X = [0 1; 1 0]
  # Z Pauli
  Z = [1 0; 0 -1]
  id = eye(2)
  return (id ⊗ X) * rho * (id ⊗ X)'
end

x = ronaldState(0.5)
ρ_ronald(p::Real) = ronald2State(p)
bennett(ρ_ronald(0))

plot(p -> fidelity(rot(ρ_ronald(p))), 0, 1, Guide.xlabel("Input F"), Guide.ylabel("Output F"), Guide.xticks(ticks=[0:0.1:1]), Guide.yticks(ticks=[0:0.1:1]))
# plot(p -> bennett(ρ_ronald(p))[1], 0, 1, Guide.xlabel("Input F"), Guide.ylabel("Output F"), Guide.xticks(ticks=[0:0.1:1]), Guide.yticks(ticks=[0:0.1:1]))
# plot(p -> fidelity(ρ_ronald(p)), 0, 1, Guide.xticks(ticks=[0:0.1:1]), Guide.yticks(ticks=[0:0.1:1]))
# typeof(DomainError  )
# println("Fidelity: ", round(fidelity(x), 3), " -> ", round(fidelity(polarize(x)), 3))
# println("Fidelity: ", round(fidelity(rot(x)), 3), " -> ", round(fidelity(polarize(rot(x))), 3))
# println("              -> ", round(fidelity(rot(polarize(x))), 3))
