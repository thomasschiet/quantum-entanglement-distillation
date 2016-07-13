print("loading packages")
tic()
print(".")
include("BEPQuantum.jl")
print(".")
using BEPQuantum
print(".")
using Iterators
println(".")
toc()

println()
print("constructing net")
print(".")
tic()

# the squares of both these nets sum up to one, so construct radial nets for both
smallNet = constructRadialNet( 2, 2)# smallNet = constructNet( 2, 3)
smallNet = constructRadialNet( 2, 2)# smallNet = constructNet( 2, 3)
print(".")
bigNet   = constructRadialNet(13, 2, 4/√6)
bigNet   = constructNet(13, 2)
print(".")
typeof(bigNet)
# the product of the small and big nets are stored in finalnet
# also add 1/6 for the last one
# 1/6 ( id_6 + c_x kron(pauli_x, id_3) + c_z kron(pauli_z, id_3) + sqrt(6) \sum_ij c_ij kron(\sigma_ij, \lambda_ij))

import Base.zero
zero(x::Type{Array{Float64, 1}}) = [0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.]

t = length(smallNet)*length(bigNet)
finalnet = Array{Float64, 1}[]
finalnet = zeros(Array{Float64, 1}, t)
i = 1
text = "0%"
println()
print(text)
for p in product(smallNet, bigNet)
  if i % 2048 == 0
    print("\b\b\b\b\b\b\b\b\b\b\b")
    text = string(round(100*i/t, 2), "%")
    print(text)
  end
  finalnet[i] = collect(tuple(1, p[1]..., p[2]...))
  # if i > 100000
    # break
  # end
  i += 1
end
println(".")
typeof(finalnet[1])

toc()

print("\n constructing basis")
tic()

# construct basis
print(".")
id2 = eye(2)
id3 = eye(3)
pauli1 = [0    1;   1   0]
pauli2 = [0 -1im; 1im   0]
pauli3 = [1    0;   0  -1]

# Gell-Mann matrices
print(".")
lambda1 =          [0    1    0;   1  0    0;   0   0  0]
lambda2 =          [0 -1im    0; 1im  0    0;   0   0  0]
lambda3 =          [1    0    0;   0 -1    0;   0   0  0]
lambda4 =          [0    0    1;   0  0    0;   1   0  0]
lambda5 =          [0    0 -1im;   0  0    0; 1im   0  0]
lambda6 =          [0    0    0;   0  0    1;   0   1  0]
lambda7 =          [0    0    0;   0  0 -1im;   0 1im  0]
lambda8 = (1/√3) * [1    0    0;   0  1    0;   0   0 -2]

# exclude all complex matrices and all matrices in the form (Gell-Mann) ⊗ id
basis = Any[]
for n in Any[id3, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8]
    y = 0
    for m in Any[id2, pauli1, pauli2, pauli3]
      r = sparse(m ⊗ n)
      if isreal(r) && ! (m == id2 && ! (n == id3))
        # cast the matrix to a real matrix
        # because it might be a complex type
        # while the complex part is 0.
        append!(basis, Any[real(r)])
      end
    end
end

# the following has to be asserted
# in order for the coefficients to
# correspond to the right basis matrix
# we're after matrices in the form
# 1/6 [id_6 + c_x (pauli_x ⊗ id_3) + c_z (pauli_z ⊗ id_3) + √6 Σ_ij c_ij (\sigma_ij ⊗ \lambda_ij)]
# This implies that, c[1] = 1/6
# c[2] ≤ 1/6
# c[3] ≤ 1/6
# c[4+] ≤ 1/√6
@assert basis[1] == id2 ⊗ id3
@assert basis[2] == pauli1 ⊗ id3
@assert basis[3] == pauli3 ⊗ id3
@assert basis[4] == pauli1 ⊗ lambda1
@assert basis[5] == pauli3 ⊗ lambda1
println(".")
toc()

# pick random
ρ = ronaldStateQutrit(1)
δ = 0.9
ϵ = 0.05

tic()
t = length(smallNet)*length(bigNet)
i = 1
text = "0%"
println("Calculating choi states")
print(text)
 for c in finalnet
  findChoiStateA(c, basis, ρ, δ, ϵ, verbose = false)
  findChoiStateB(c, basis, ρ, δ, ϵ, verbose = false)
  if i % 4 == 0
    print("\b\b\b\b\b\b\b\b\b\b\b")
    text = string(round(100*i/t, 2), "%")
    print(text)
  end
  i += 1
end
toc()

getRandomC(cs) = cs[round(Int, rand() * length(finalnet))]
function findRandomChoiB(cs, basis, ρ, δ, ϵ)
  c = getRandomC(cs)
  println(round(c, 3))
  return findChoiStateB(c, basis, ρ, δ, ϵ)
end
function findRandomChoiA(cs, basis, ρ, δ, ϵ)
  c = getRandomC(cs)
  println(round(c, 3))
  return findChoiStateA(c, basis, ρ, δ, ϵ)
end

@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time findRandomChoiB(finalnet, basis, ρ, δ, ϵ)
@time (problem, C_A, C_1A, C_0A) = findRandomChoiA(finalnet, basis, ρ, δ, ϵ)
@time (problem, C_B, C_1B, C_0B) = findRandomChoiB(finalnet, basis, ρ, δ, ϵ)

println(fidelityOfChoistate(C_1A, C_1B, ρ))
println(fidelityOfChoistate(C_1B, C_1A, ρ))
println(fidelityOfChoistate(C_1A, C_1A, ρ))
println(fidelityOfChoistate(C_1B, C_1B, ρ))

fidelityOfChoistate(C_1A::Convex.Variable, C_1B::Convex.Variable, ρ::Union{Matrix, SparseMatrixCSC}) = fidelityOfChoistate(C_1A.value, C_1B.value, ρ)
function fidelityOfChoistate(C_1A::Union{Matrix, SparseMatrixCSC}, C_1B::Union{Matrix, SparseMatrixCSC}, ρ::Union{Matrix, SparseMatrixCSC})
  e0 = eVec(2, 1)
  e1 = eVec(2, 2)

  epr = (e0 ⊗ e0 + e1 ⊗ e1)
  epr = epr * epr'

  return trace((epr ⊗ ρ') * (C_1A⊗C_1B))
end
