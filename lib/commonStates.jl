#
# Defines some commonly used states in quantum information
#
export wernerState, ronaldState, ronald2State, epr

include("wernerState.jl")
include("ronaldState.jl")
include("ronald2State.jl")
include("ronaldStateQutrit.jl")
include("eVec.jl")

epr = (kron(eVec(2,1),eVec(2,1)) + kron(eVec(2,2),eVec(2,2)))/sqrt(2);
singlet = (kron(eVec(2,1),eVec(2,2)) - kron(eVec(2,2),eVec(2,1)))/sqrt(2);
ρEpr = epr * epr';

bell = zeros(4,4);
bell[1,1:4] = epr;
bell[2,1:4] = (kron(eVec(2,1),eVec(2,1)) -  kron(eVec(2,2),eVec(2,2)))/sqrt(2);
bell[3,1:4] = (kron(eVec(2,1),eVec(2,2)) +  kron(eVec(2,2),eVec(2,1)))/sqrt(2);
bell[4,1:4] = singlet;

bS1 = epr*epr';
bS2 = bell[2,1:4]'* bell[2,1:4];
bS3 = bell[3,1:4]'* bell[3,1:4];
bS4 = bell[4,1:4]'*bell[4,1:4];

bellS = [bS1; bS2; bS3; bS4];

PS = bS1 + bS2 + bS3;
PA = bS4;


v0 = eVec(2,1);
v1 = eVec(2,2);

v00 = kron(v0,v0);
v11 = kron(v1,v1);
v01 = kron(v0,v1);
v10 = kron(v1,v0);
