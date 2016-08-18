using Convex
using SCS

export PPTprogrammeNoTwirling1Ext
export PPTprogrammeNoTwirling2Ext
export PPTprogrammeNoTwirling1ExtPermSym
export PPTprogrammeNoTwirling1ExtPermSymOnlySym

"""
 Implements the Rains bound for distillable entanglement
 with an allowed failure probability (which in turn may give higher fidelity)
 and with 1-extension

 Inputs:
 rho quantum state
 nA  dimension Alice
 nB  dimension Bob (assumed to be the same for now)
 K   output dimension
 δ   maximum allowed failure probability
 n   number of input copies

 Outputs:
 problem problem object given by Convex. Optimal value can be obtained
         by problem.optval
"""
function PPTprogrammeNoTwirling1Ext(rho, nA, nB, n, K, δ; verbose = true)

  # Check whether rho is a quantum state
  @assert isQuantumState(rho)

  # Check whether dimensions match
  (d, db) = size(rho)
  @assert d == nA*nB

  # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  if n > 1
    rhoSorted = sortAB(rho, 2, n)
  else
    rhoSorted = rho
  end

  # define the variables W
  W_AB1B2 = Semidefinite(2^(3 * (n + 1)))

  # Define the EPR pair:
  epr = (eVec(2, 1) ⊗ eVec(2, 1) + eVec(2, 2) ⊗ eVec(2, 2))/√2
  epr = epr * epr'

  # dimensions of W
  dims = [2, 2^n, # Ahat A'
          2, 2^n, # B_1hat B_1'
          2, 2^n] # B_2hat B_2'
  # Choi with first B system
  W_AB1 = partialtrace(W_AB1B2, [5, 6], dims)

  # Choi with second B system
  W_AB2 = partialtrace(W_AB1B2, [3, 4], dims)

  # define the objective
  problem = maximize(nA * nB * trace(swap(epr ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1))

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_AB1 , [1, 3], [2, 2^n, 2, 2^n]) ⪰ 0

    # constrain probability of success
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≤ δ
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≥ 0

    ptranspose(W_AB1) ⪰ 0

    # Choi should be equal if we trace out one of the extensions
    W_AB1 == W_AB2
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose, eps = 1e-3))

  # Output
  Psuccess = nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(W_AB1B2.value, [5, 6], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess)
end

function PPTprogrammeNoTwirling2Ext(rho, nA, nB, n, K, δ; verbose = true)

  # Check whether rho is a quantum state
  @assert isQuantumState(rho)

  # Check whether dimensions match
  (d, db) = size(rho)
  @assert d == nA*nB

  # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  if n > 1
    rhoSorted = sortAB(rho, 2, n)
  else
    rhoSorted = rho
  end

  # define the variables W
  W_AB1B2B3 = Semidefinite(2^(4 * (n + 1)))

  # Define the EPR pair:
  epr = (eVec(2, 1) ⊗ eVec(2, 1) + eVec(2, 2) ⊗ eVec(2, 2))/√2
  epr = epr * epr'

  # dimensions of W
  dims = [2, 2^n, # Ahat A'
          2, 2^n, # B_1hat B_1'
          2, 2^n, # B_2hat B_2'
          2, 2^n] # B_3hat B_3'
  # Choi with first B system
  W_AB1 = partialtrace(W_AB1B2B3, [5, 6, 7, 8], dims)

  # Choi with second B system
  W_AB2 = partialtrace(W_AB1B2B3, [3, 4, 7, 8], dims)

  W_AB3 = partialtrace(W_AB1B2B3, [5, 6, 7, 8], dims)

  # define the objective
  problem = maximize(nA * nB * trace(swap(epr ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1))

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_AB1 , [1, 3], [2, 2^n, 2, 2^n]) ⪰ 0

    # constrain probability of success
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≤ δ
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≥ 0

    ptranspose(W_AB1) ⪰ 0

    # Choi should be equal if we trace out one of the extensions
    W_AB1 == W_AB2
    W_AB2 == W_AB3
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose))

  # Output
  Psuccess = nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(W_AB1B2B3.value, [5, 6, 7, 8], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess)
end

function PPTprogrammeNoTwirling1ExtPermSymOnlySym(rho, nA, nB, n, K, δ; verbose = true)

  # Check whether rho is a quantum state
  @assert isQuantumState(rho)

  # Check whether dimensions match
  (d, db) = size(rho)
  @assert d == nA*nB

  # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  if n > 1
    rhoSorted = sortAB(rho, 2, n)
  else
    rhoSorted = rho
  end

  # define the variables W
  Ws = Semidefinite(288)
  Wa = spzeros(224, 224)

  Pcb = getPcb()

  # Define the EPR pair:
  epr = (eVec(2, 1) ⊗ eVec(2, 1) + eVec(2, 2) ⊗ eVec(2, 2))/√2
  epr = epr * epr'

  # dimensions of W
  dims = [2, 2^n, # Ahat A'
          2, 2^n, # B_1hat B_1'
          2, 2^n] # B_2hat B_2'

  W_AB1B2 = Pcb * (Ws ⊕ Wa) * Pcb'

  # Choi with first B system
  W_AB1 = partialtrace(W_AB1B2, [5, 6], dims)

  # Choi with second B system
  W_AB2 = partialtrace(W_AB1B2, [3, 4], dims)

  # define the objective
  # problem = maximize(nA * nB * trace(swap(epr ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1))
  X =  Pcb' * nA * nB *(swap(epr ⊗ rhoSorted', [2, 3], dim = [2, 2, 2^n, 2^n]) ⊗ eye(2 * 2^n)) * Pcb
  problem = maximize(trace(X[1:288, 1:288] * Ws))

  # 5/0.6 = 4.12 / p
  #

  Y = Pcb' * nA * nB * (swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) ⊗ eye(2 * 2^n)) * Pcb
  p_succ = trace(Y[1:288, 1:288] * Ws)

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_AB1 , [1, 3], [2, 2^n, 2, 2^n]) ⪰ 0

    # constrain probability of success
    p_succ ≤ δ
    p_succ ≥ 0

    ptranspose(W_AB1) ⪰ 0
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose, eps = 1e-4))

  # Output
  println(trace(Y[1:288, 1:288] * Ws.value))
  Psuccess = nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(Pcb * (Ws.value ⊕ Wa) * Pcb', [5, 6], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess)
end

function PPTprogrammeNoTwirling1ExtPermSym(rho, nA, nB, n, K, δ; verbose = true)

  # Check whether rho is a quantum state
  @assert isQuantumState(rho)

  # Check whether dimensions match
  (d, db) = size(rho)
  @assert d == nA*nB

  # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  if n > 1
    rhoSorted = sortAB(rho, 2, n)
  else
    rhoSorted = rho
  end

  # define the variables W
  Ws = Semidefinite(288)
  Wa = Semidefinite(224, 224)

  Pcb = getPcb()

  # Define the EPR pair:
  epr = (eVec(2, 1) ⊗ eVec(2, 1) + eVec(2, 2) ⊗ eVec(2, 2))/√2
  epr = epr * epr'

  # dimensions of W
  dims = [2, 2^n, # Ahat A'
          2, 2^n, # B_1hat B_1'
          2, 2^n] # B_2hat B_2'

  W_AB1B2 = Pcb * (Ws ⊕ Wa) * Pcb'

  # Choi with first B system
  W_AB1 = partialtrace(W_AB1B2, [5, 6], dims)

  # Choi with second B system
  W_AB2 = partialtrace(W_AB1B2, [3, 4], dims)

  # define the objective
  problem = maximize(nA * nB * trace(swap(epr ⊗ rhoSorted', [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1))

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_AB1 , [1, 3], [2, 2^n, 2, 2^n]) ⪰ 0

    # constrain probability of success
    nA * nB * trace(swap(eye(4) ⊗ rhoSorted', [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≤ δ
    nA * nB * trace(swap(eye(4) ⊗ rhoSorted', [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≥ 0

    ptranspose(W_AB1) ⪰ 0

    # Choi should be equal if we trace out one of the extensions
    # W_AB1 == W_AB2
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose, eps = 1e-3))

  # Output
  Psuccess = nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(Pcb * (Ws.value ⊕ Wa.value) * Pcb', [5, 6], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess, Wa)
end

function getPcb()
  v0 = spzeros(8, 1)
  v = Any[]
  for i = 1:8
    x = copy(v0)
    x[i] = 1
    push!(v, x)
  end

  u = Any[]
  for i = 1:8
    x = kron( v[i] , v[i] )
    push!(u, x)
  end

  # symmetric
  for j = 1:7
    for i = 1:8-j
      x = ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) + kron( v[i+j] , v[j] ) )
      push!(u, x)
    end
  end

  # asymmetric
  for j = 1: 7
    for i = 1:8-j
      x = ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) - kron( v[i+j] , v[j] ))
      push!(u, x)
    end
  end

  Pcb = spzeros(512, 512)
  for i = 1:64
    for j = 1:8
      Pcb[:, 64 * (j - 1) + i] = kron(v[j], u[i])
    end
  end

  return Pcb
end
