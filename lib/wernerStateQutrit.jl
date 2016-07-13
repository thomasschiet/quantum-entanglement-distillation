module Quantum

function wernerStateQutrit(p)
  e0 = [1 0 0]'
  e1 = [0 1 0]'
  e2 = [0 0 1]'


  epr = 1/√3 * (e0 ⊗ e0 + e1 ⊗ e1 + e2 ⊗ e2)
  id = eye(9)

  return epr*epr' + id/9
end

wernerStateQutrit(0.2)

end
