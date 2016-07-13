include("commonStates.jl")

function bellDiag(state)
  result = zeros(4)*zeros(4)'
  for i = 1:4
    for j = 1:4
      result[i, j] = (bell[i,1:4] * state * bell[j,1:4]')[1];
    end
  end
  return result
end

function deutsch(state)
  return bellDiag(state)
end

ρ = wernerState(p)
result[1, 1] = (bell[1,1:4] * ρ * bell[1, 1:4]')[1]

result = zeros(4)*zeros(4)'


p = 1
ρ = ronald2State(p)
# print(deutsch(ρ))
# eprFidelity(ρ)
# print(ρ)
