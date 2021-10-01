% This function is an n-dimensional generalization of mySpread. In this
% function there are no time steps; just one interaction for each
% individual.
function f = diseasespread(f0, dt)

global mu beta

f = f0;

% Population, susceptibles and infectious in each node
nodePop = sum(f0);
nS = f0(1,:);
nI = f0(2,:);

% Get indices of nodes where the population is not zero. These are the
% indices that are treated in this function.
indicies = find(nodePop ~= 0);
% This is the maximum total population of any node. This number is needed
% for knowing how long vectors holding the individuals in a node must be.
maxPop = max(nodePop);

% Create an array where col index denotes a specific node and row index
% denotes a particular individual
S = zeros(maxPop, length(indicies));
I = zeros(maxPop, length(indicies));

r = rand(maxPop, length(indicies));

for i = 1:length(indicies)
    S(1 : nS(i), i) = 1;
    I((nS(i) + 1 : nodePop(i)), i) = 1;
end
pItoS = mu*dt;
pStoI = repmat(1-(1-beta*dt).^nI(indicies), maxPop, 1);
%pStoI = repmat(pStoI, length(indicies), 1)


ItoS = (I.*r<pItoS).*(I.*r>0);
StoI = (S.*r<pStoI).*(S.*r>0);

S = S - StoI + ItoS;
I = I - ItoS + StoI;

f(1,indicies) = sum(S);
f(2,indicies) = sum(I);
% Infectious and susceptible vectors

% S(1:f0(1)) = 1
% I = zeros(1,sum(f0));
% I((f0(1)+1):sum(f0)) = 1;
