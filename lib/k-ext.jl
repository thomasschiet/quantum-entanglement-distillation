using Convex
using SCS

export PPTprogrammeNoTwirling1Ext

"""
 Implements the Rains bound for distillable entanglement
 with an allowed failure probability (which in turn may give higher fidelity)
 and with 1-extension

 Inputs:
 rho quantum state
 nA  dimension Alice
 nB  dimension Bob (assumed to be the same for now)
 K   output dimension
 delta   maximum allowed failure probability
 n   number of input copies

 Outputs:
 problem problem object given by Convex. Optimal value can be obtained
         by problem.optval
"""
function PPTprogrammeNoTwirling1Ext(rho, nA, nB, n, K, delta verbose = true)

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

  # define the objective
  problem = maximize(nA * nB * trace(swap(epr ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(W, [5, 6], dims)))

  # Choi with first B system
  W_AB1 = partialtrace(W, [5, 6], dims)

  # Choi with second B system
  W_AB2 = partialtrace(W, [3, 4], dims)

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_AB1 , [1, 3], [2, 2^n, 2, 2^n]) ⪰ 0

    # constrain probability of success
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≤ delta
    nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * W_AB1) ≥ 0

    partialtranspose(W_AB1) ⪰ 0

    # Choi should be equal if we trace out one of the extensions
    W_AB1 == W_AB2
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose))

  # Output
  Psuccess = nA * nB * trace(swap(eye(4) ⊗ transpose(rhoSorted), [2, 3], dim = [2, 2, 2^n, 2^n]) * partialtrace(W.value, [5, 6], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess)
end
