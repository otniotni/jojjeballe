close all
%% Task 1
close all
clear all
global mu beta f0
% Set time
tmax = 15000;
tspan = [0 tmax];

% Set initial condition
I0 = 35;
S0 = 1000;
% I0 = 5;
% S0 = 95;
f0 = [S0;I0];
% Set parameters
% Note: For mu > Ntot*beta/2, the number of infected individuals reach a steady
%       state where S > I
mu = 1e-3; % Rate of an individual recovering from the infection
beta = 1.5*mu/sum(f0); % Rate of an individual becoming infected due to proximity 
          % to another infected individual

% Fixed point
Sfix = mu/beta;
Ifix = sum(f0)-Sfix;


% Set ode-related things
rtol = 1e-6;
atol = 1e-4;
opts = odeset('RelTol',rtol,'AbsTol',atol);

odeII = @(t,f) [mu*f(2)-beta*f(1).*f(2);-mu*f(2)+beta*f(1).*f(2)];

% Solve ode
solII = ode45(odeII,tspan,f0,opts);

% Calculate relative populations
Srel = solII.y(1,:)./(solII.y(1,:)+solII.y(2,:));
Irel = solII.y(2,:)./(solII.y(1,:)+solII.y(2,:));

% Plot(not relative)
figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
plot(solII.x,solII.y)
yline(Sfix)
yline(Ifix)
legend('S(t)','I(t)','S^{*}','I^{*}')


% Task 2

% Good way to find dt?
dt = 5e-2;
t = tspan(1):dt:tspan(end);
f2 = zeros(2, length(t));
f2(:,1) = f0;
tic;[t,f] = mySpread(tspan, dt);toc
%f = diseasespread(f0, dt)
figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
plot(t,f)
% plot(t,f2)
plot(solII.x,solII.y)
yline(Sfix)
yline(Ifix)
legend('S_{stochastic1}','I_{stochastic1}',...
    'S_{stochastic2}','I_{stochastic2}','S','I')
% Plot(relative)
figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
plot(solII.x,Srel)
plot(solII.x,Irel)
plot(t, f/sum(f0))
yline(Sfix/sum(f0),'-','$\rho_S = \frac{\mu}{\beta(I_0+S_0)}$','Interpreter','latex','LineWidth',4,'FontSize',28,'LabelHorizontalAlignment','left')
yline(Ifix/sum(f0),'-','$\rho_I = 1 - \frac{\mu}{\beta(I_0+S_0)}$','Interpreter',...
    'latex','LineWidth',4,'FontSize',28,'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left')
legend('$\rho_{S}$','$\rho_{I}$','$\rho_{S, stochastic}$','$\rho_{I,stochastic}$','Interpreter','latex')
xlabel('Time')
ylabel('Relative population')



%% Task 4
% clear all
close all
global neighbourCount neighbours n
% Load data
pownet = load('pownet_N100_gamma-2.5.mat');
W = pownet.W;
W(W ~= 0) = 1;
pung = graph(W);

% Number of nodes
n = size(W,1);
% Randomly distributed populations
rho = 50; % Average population per node
var = floor(rho/5); % Spread
% Population distribution: This becomes increasingly fucked the larger var
% is. Should maybe find a better way to do this in order to preserver the
% mean while still having a deviation.
% pop = var*randi([-rho rho],1,100) + rho;
% mean1 = mean(pop)
pop = randomPop(rho,n);

% Degree stuff
k = sum(W); % Degree in node
kmax = max(k); % Maximum degree in the system
kmean = mean(sum(W)); % Mean node degree

for ne = 1:n
    nbs = find(W(ne,:));
    neighbourCount(ne) = length(nbs);
    neighbours(ne,1:neighbourCount(ne)) = nbs;
end


tmax = 1000;

% Allocate population mat
popMat = zeros(tmax,length(pop));

for i = 1:tmax
    pop = mySpread2mod(W, pop); % Calculate new population
    popMat(i,:) = pop; % Insert latest population into popMat
end
% Mean of the populations for each node
popMean = mean(popMat,1);

% Fit popMean against k
p = polyfit(k, popMean, 1);

% Create fitted function and theoretical function and calculate difference
linFit = @(x) p(1)*x + p(2);
realFun = @(x) x*rho/kmean;
difference = abs(p(1)-rho/kmean)/(rho/kmean);

figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
scatter(k, popMean)
fplot(linFit)
fplot(realFun)
ylabel('Average number of individuals')
xlabel('Node degree')
axis([0 kmax 0 max(popMean)])
legend('Experimental data','Linear fit','Equation(3)')


%% Task 4 continued(for different gamma)
close all
pownet25 = load('pownet_N100_gamma-2.5.mat');
pownet30 = load('pownet_N100_gamma-3.0.mat');
pownet35 = load('pownet_N100_gamma-3.5.mat');
pownet40 = load('pownet_N100_gamma-4.0.mat');
W25 = pownet25.W;
W25(W25 ~= 0) = 1;
W30 = pownet30.W;
W30(W30 ~= 0) = 1;
W35 = pownet35.W;
W35(W35 ~= 0) = 1;
W40 = pownet40.W;
W40(W40 ~= 0) = 1;

% Number of nodes
n = size(W,1);
% Randomly distributed populations
rho = 50; % Average population per node
var = floor(rho/5); % Spread

pop25 = randomPop(rho,n);
pop30 = randomPop(rho,n);
pop35 = randomPop(rho,n);
pop40 = randomPop(rho,n);

k25 = sum(W25); % Degree in node
k30 = sum(W30); % Degree in node
k35 = sum(W35); % Degree in node
k40 = sum(W40); % Degree in node
kmax25 = max(k25); % Maximum degree in the system
kmax30 = max(k30); % Maximum degree in the system
kmax35 = max(k35); % Maximum degree in the system
kmax40 = max(k40); % Maximum degree in the system
kmean25 = mean(sum(W25)); % Mean node degree
kmean30 = mean(sum(W30)); % Mean node degree
kmean35 = mean(sum(W35)); % Mean node degree
kmean40 = mean(sum(W40)); % Mean node degree

tmax = 1000;

% Allocate population mat
popMat25 = zeros(tmax,length(pop25));
popMat30 = zeros(tmax,length(pop30));
popMat35 = zeros(tmax,length(pop35));
popMat40 = zeros(tmax,length(pop40));

for i = 1:tmax
    pop25 = mySpread2(W25, pop25); % Calculate new population
    pop30 = mySpread2(W30, pop30); % Calculate new population
    pop35 = mySpread2(W35, pop35); % Calculate new population
    pop40 = mySpread2(W40, pop40); % Calculate new population
    popMat25(i,:) = pop25; % Insert latest population into popMat
    popMat30(i,:) = pop30; % Insert latest population into popMat
    popMat35(i,:) = pop35; % Insert latest population into popMat
    popMat40(i,:) = pop40; % Insert latest population into popMat
end
popMean25 = popMat25(end, :);
popMean30 = popMat30(end, :);
popMean35 = popMat35(end, :);
popMean40 = popMat40(end, :);
% popMean25 = mean(popMat25(tmax-1:tmax,:),1); % Mean of the populations for each node
% popMean30 = mean(popMat30(tmax-1:tmax,:),1); % Mean of the populations for each node
% popMean35 = mean(popMat35(tmax-1:tmax,:),1); % Mean of the populations for each node
% popMean40 = mean(popMat40(tmax-1:tmax,:),1); % Mean of the populations for each node
p25 = polyfit(k25, popMean25, 1);
p30 = polyfit(k30, popMean30, 1);
p35 = polyfit(k35, popMean35, 1);
p40 = polyfit(k40, popMean40, 1);
% Create fitted function and theoretical function and calculate difference
linFit25 = @(x) p25(1)*x + p25(2);
linFit30 = @(x) p30(1)*x + p30(2);
linFit35 = @(x) p35(1)*x + p35(2);
linFit40 = @(x) p40(1)*x + p40(2);
realFun25 = @(x) x*rho/kmean25;
realFun30 = @(x) x*rho/kmean30;
realFun35 = @(x) x*rho/kmean35;
realFun40 = @(x) x*rho/kmean40;

difference25 = abs(p25(1)-rho/kmean25)/(rho/kmean25)
difference30 = abs(p30(1)-rho/kmean30)/(rho/kmean30)
difference35 = abs(p35(1)-rho/kmean35)/(rho/kmean35)
difference40 = abs(p40(1)-rho/kmean40)/(rho/kmean40)
figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
% scatter(k, popMean)
sz = 150;
fplot(realFun25,'LineWidth',3)
fplot(realFun30,'LineWidth',3)
fplot(realFun35,'LineWidth',3)
fplot(realFun40,'LineWidth',3)
fplot(linFit25,'--','LineWidth',3)
fplot(linFit30,'--','LineWidth',3)
fplot(linFit35,'--','LineWidth',3)
fplot(linFit40,'--','LineWidth',3)
% scatter(k25, popMean25,sz, 'o','filled','MarkerEdgeColor','black')
% scatter(k30, popMean30,2*sz, '<','filled','MarkerEdgeColor','black')
% scatter(k35, popMean35,sz, 's','filled','MarkerEdgeColor','black')
% scatter(k40, popMean40,sz, 'd','filled','MarkerEdgeColor','black')
ylabel('Average number of individuals')
xlabel('Node degree')
axis([0 kmax 0 max(popMean)])
legend('$\gamma = -2.5: \rho_k = \frac{\rho}{\left\langle k \right\rangle}k$',...
    '$\gamma = -3.0: \rho_k = \frac{\rho}{\left\langle k \right\rangle}k$',...
    '$\gamma = -3.5: \rho_k = \frac{\rho}{\left\langle k \right\rangle}k$',...
    '$\gamma = -4.0: \rho_k = \frac{\rho}{\left\langle k \right\rangle}k$',...
    '$\gamma = -2.5: Stochastic$', '$\gamma = -3.0: Stochastic$',...
    '$\gamma = -3.5: Stochastic$','$\gamma = -4.0: Stochastic$','Interpreter','latex')
    


% \rho_k = k*(\rho/kmean) slope should be rho/kmean

%% Task 5a
global betafunc status
filename = '';
status = 1;
beta = 1; mu = 2;
betafunc = 1;
Nspan = 2:4;

gammaspan = 2.5:2.5;
rhospan = .5:.1:2;
iterations = 1:10;
% avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations, 1, 1, filename);

%% Task 5a cont
close all
plot56(load('averageIrelativev1.mat'), rhospan);
%% Task 5b
rhospan = 1.5:.1:3;
filename = 'noSus';
status = 2;
% avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations, 2, 1, filename);
%% Task 5b cont
close all
plot56(load('averageIrelative-noSusv1.mat'), rhospan)
%% Task 6
filename = 'newBeta';
beta = 1; mu = .5;
betafunc = @(rhoi) 1/rhoi;
status = 3;
rhospan = .5:.1:2;
iterations = 1:10;
% avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations, 3, 0, filename);
%% Task 6 cont
plot56(load('averageIrelative-newBetav1.mat'), rhospan)
%% Task 6 bonus
% a) We try to reduce the spreading by not letting infected individuals
% travel.
filename = 'noInf';
% betafunc = 1;
status = 4;
avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations, status, 1, filename);
% b) We try to reduce the spreading by completely shutting down travel to
% and from nodes with an Irelative above a given threshold.
filename = 'quarantine';
rhospan = .5:.1:2;
status = 5;
% avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations, status, 1, filename);
%% Task 6 bonus cont
plot56(load('averageIrelative-quarantinev1.mat'),rhospan)
%% Task 6 bonus cont
plot56(load('averageIrelative-noInfv1.mat'),rhospan)