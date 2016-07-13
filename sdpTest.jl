# Test an SDP

using Convex

# dimension Alice and Bob
nA = 2;
nB = 2;

n = nA * nB;
id = eye(n,n);

# define the choi state variable
choiState = Variable(n,n);

A = randn(n,n);

problem = maximize(trace(A * choiState));

problem.constraints += ([choiState in :SDP]);
problem.constraints += ([trace(choiState) == 1]);


solve!(problem)

problem.status
problem.optval

