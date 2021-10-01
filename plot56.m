function [h1,h2,h3] = plot56(Idata, rhospan)

Idata = Idata.avgIrel;
I100rel = Idata(:,:,1);
I1000rel = Idata(:,:,2);
I10000rel = Idata(:,:,3);
% This plot shows a scatter plot of the last 200 timesteps of Irelative for
% rho=0.5:0.1:2 for N = {100,1000,10000}.
figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
h1 = scatter(rhospan, I100rel,100,'red', 'd','filled');
h2 = scatter(rhospan, I1000rel, 100,'blue', 'p', 'filled');
h3 = scatter(rhospan, I10000rel,100,'black', 'filled');
xlabel('$\rho$','Interpreter','latex')
ylabel('$\frac{\rho_I}{\rho}$','Interpreter','latex','Rotation',0,...
    'HorizontalAlignment','right')
legend([h1(1),h2(1),h3(1)],'$N=10^{2}$','$N=10^{3}$'...
    ,'$N=10^{4}$','Interpreter','latex')
axis([rhospan(1) rhospan(end) 0 0.6])