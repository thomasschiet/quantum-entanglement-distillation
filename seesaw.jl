include("quantum.jl")

using Convex
using SCS

function seesawBfixed(ρ_input, nCopies, nA, nB, K, δ_min, δ_max, C_B)

  epr = (kron(eVec(2,1),eVec(2,1)) + kron(eVec(2,2),eVec(2,2)))/sqrt(2);
  max = epr*epr';


  # make copies of the input state
  ρ_input = sortAB(copies(ρ_input, nCopies), nCopies, nCopies)
  ρ_input_transposed = complexToReal(ρ_input')

  if (!(1 > trace(C_B) > 0))
    error("C_B has invalid trace: ", trace(C_B))
    return
  end

  (d,db) = size(ρ_input);
  if (d != nA*nB)
    error("Input state doesn't match given dimensions." , d , "≠", nA*nB);
    return;
  end

  id = complexToReal(eye(d,d));

  C_A = Variable(2d, 2d)
  C_B = complexToReal(C_B)

  problem = maximize(nA * nB * trace(kron(epr, ρ_input_transposed) * kron(C_A, C_B)) / 2);

  # constrain the probability of success
  problem.constraints += ([nA * nB * trace(ρ_input_transposed * kron(C_A, C_B)) ≥ δ_min])
  problem.constraints += ([nA * nB * trace(ρ_input_transposed * kron(C_A, C_B)) ≤ δ_max])

  # make sure that C_A is a valid Choi state
  problem.constraints += ([C_A in :SDP]);
  problem.constraints += ([trace(C_A) ≥ 0]);
  problem.constraints += ([id/nA - C_A in :SDP]);

  solve!(problem, SCSSolver(verbose=true));
end

p = 0.7
ρ_input = wernerState(p)
nCopies = 2
nA = 4
nB = 4
K = 2
δ_min = 0
δ_max = 0.8
C_B = kron(eye(4)/4.1, eye(4)/4.1)
print(seesawBfixed(ρ_input, nCopies, nA, nB, K, δ_min, δ_max, C_B))
