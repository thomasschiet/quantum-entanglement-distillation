using Convex
using SCS

TOL = 1e-3

facts("Permute systems") do
  # check examples from http://www.qetlab.com/PermuteSystems
  X = reshape(collect(1:16), (4, 4))'
  @fact PermuteSystems(X, [2; 1]) --> [1 3 2 4; 9 11 10 12; 5 7 6 8; 13 15 14 16]
  @fact PermuteSystems(X, [2; 1], dim = [2; 2]) --> [1 3 2 4; 9 11 10 12; 5 7 6 8; 13 15 14 16]

  @fact PermuteSystems(X, [1; 2]) --> X
  @fact PermuteSystems(X, [1; 2], dim = [2; 2]) --> X

  X = rand(4, 4)
  @fact PermuteSystems(PermuteSystems(X, [2; 1]), [2; 1]) --> X

  A = rand(4, 4)
  B = rand(4, 4)
  @fact PermuteSystems(A ⊗ B, [2; 1]) --> B ⊗ A

  A = rand(4, 4)
  B = rand(8, 8)
  @fact PermuteSystems(A ⊗ B, [2; 1], dim = [4; 8]) --> B ⊗ A

  X = reshape(collect(1:64), (8, 8))'
  Y = [     1     5     2     6     3     7     4     8;
    33    37    34    38    35    39    36    40;
     9    13    10    14    11    15    12    16;
    41    45    42    46    43    47    44    48;
    17    21    18    22    19    23    20    24;
    49    53    50    54    51    55    52    56;
    25    29    26    30    27    31    28    32;
    57    61    58    62    59    63    60    64;]
  @fact PermuteSystems(X, [2; 3; 1]) --> Y

  # Test if permute systems works in Convex.jl.
  # This is an revised version of the following test
  # that was taken from the Convex.jl source.
  # https://github.com/JuliaOpt/Convex.jl/blob/master/test/test_affine.jl
  #
  # X = Variable(2, 2)
  # c = ones(2, 1)
  # p = minimize(c' * X' * c, [X >= ones(2, 2)])
  # @fact vexity(p) --> AffineVexity()
  # solve!(p)
  # @fact p.optval --> roughly(4, TOL)
  # @fact evaluate(c' * X' * c)[1] --> roughly(4, TOL)
  # X = Variable(4, 4)
  # c = ones(4, 1)
  # p = minimize(c' * PermuteSystems(X, [2; 1]) * c, [X >= ones(4, 4)])
  # @fact vexity(p) --> AffineVexity()
  # solve!(p, SCSSolver(verbose = false))
  # @fact p.optval --> roughly(16, TOL)
  # @fact evaluate(c' * PermuteSystems(X, [2; 1]) * c)[1] --> roughly(16, TOL)
end
