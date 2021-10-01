function pops = myCombinedSpread(dt, pops, W, status)

global n f0 betafunc beta
betasave = beta;
tspan = [0 dt];
% (a) Spread disease/get susceptible 
tstart = tic;
for i = 1:n
    f0 = pops(:,i);
    if(isa(betafunc,'function_handle'))
        if(sum(f0))
            beta = betafunc(sum(f0))*beta;
        else
            beta = 1;
        end
    end
    [t,f] = mySpread(tspan, dt);
    pops(:,i) = f(:,length(t));
    beta = betasave;
end
tend1 = toc(tstart);
tstart = tic;
% (b) Diffuse along the network
if(status == 1 || status == 3 || status == 5)
    pops = mySpread2mod(W, pops);
elseif(status == 2)
    % Susceptibles not allowed to diffuse
    pops(2,:) = mySpread2mod(W, pops(2,:));
elseif(status == 4)
    % Infectious not allowed to diffuse
    pops(1,:) = mySpread2mod(W, pops(1,:));
end
tend2 = toc(tstart);
algorithm1 = tend1;
algorithm2 = tend2;
