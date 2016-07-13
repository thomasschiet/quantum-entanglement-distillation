export bellState

function bellState(i)
  epr = (kron(eVec(2,1),eVec(2,1)) + kron(eVec(2,2),eVec(2,2)))/sqrt(2);
  singlet = (kron(eVec(2,1),eVec(2,2)) - kron(eVec(2,2),eVec(2,1)))/sqrt(2);
  ρEpr = epr * epr';

  bell = zeros(4,4);
  bell[1,1:4] = epr;
  bell[2,1:4] = (kron(eVec(2,1),eVec(2,1)) -  kron(eVec(2,2),eVec(2,2)))/sqrt(2);
  bell[3,1:4] = (kron(eVec(2,1),eVec(2,2)) +  kron(eVec(2,2),eVec(2,1)))/sqrt(2);
  bell[4,1:4] = singlet;

  if i == 1
    return epr*epr';
  elseif i == 2
    return bell[2,1:4]'* bell[2,1:4];
  elseif i == 3
    return bell[3,1:4]'* bell[3,1:4];
  elseif i == 4
    return bell[4,1:4]'*bell[4,1:4];
  end
end
