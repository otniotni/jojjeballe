% Function [t,f] = mySpread()
% ---
% Input:        tspan: vector containing first and last times
%                  dt: time step
%
% Output:
%                   t: vector of times
%                   f: vector of function data

function [t,f] = mySpread(tspan, dt)

global mu beta f0

% Time vector
t = tspan(1):dt:tspan(end);

f = zeros(2,length(t));
f(:,1) = f0;
if(sum(f0) == 0) 
    return
end

% Infected and susceptible vectors
S = zeros(1,sum(f0));
S(1:f0(1)) = 1;
I = zeros(1,sum(f0));
I((f0(1)+1):sum(f0)) = 1;


for i = 2:length(t)
    kIS = mu*dt; 
    kSI = 1-(1-beta*dt)^f(2,i-1);
    ktot = kIS+kSI;
    pItoS = kIS;
    pStoI = kSI;
    pItoS + pStoI;
    
    % Random numbers
    r = rand(1, sum(f0));
    % Indices of individuals that should go from infected/susceptible
    % to susceptible/infected, respectively
    ItoS = (I.*r<pItoS).*(I.*r>0);
    StoI = (S.*r<pStoI).*(S.*r>0);
    
    % Updating susceptible and infected individuals
    S = S - StoI + ItoS;
    I = I - ItoS + StoI;
    
    % Saving number of susceptible and infected individuals in f
    f(1,i) = sum(S);
    f(2,i) = sum(I);
end



