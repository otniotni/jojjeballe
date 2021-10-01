function avgIrel = diseasesim(Nspan, gammaspan, rhospan, iterations,...
    status, s, filename)
close all
global neighbours neighbourCount n
tmax = 500;
dt = 5e-2;
span = (tmax-100):tmax;
Irelative = zeros(1, tmax);
tol = 5e-2;


figure('DefaultAxesFontSize',24,'DefaultLineLineWidth',2)
hold on
axis([0 tmax 0 1])
% avgIrel is the mean I/rho over the last 100 seconds of the simulation,
% provided that I hasn't gone extinct. Note that it is a 3D matrix of order
% ({iterations},{rho},{N}).
avgIrel = zeros(length(iterations), length(rhospan), length(Nspan));
nidx = 0;
for pow = Nspan
    nidx = nidx + 1;
    N = 10^pow;
    for gamma = gammaspan
        % Read the network of order N as gathered from a power law
        % distribution of order gamma
        mynet = sprintf('pownet_N%d_gamma-%.1f.mat', N, gamma);
        pownet = load(mynet);
        W = pownet.W;
        W(W ~= 0) = 1;
        n = size(W,1);
        % Constructing the matrix neighbours which has all the neighbours'
        % indicies as column vectors and the vector neighbourCount which
        % gives the number of neighbours for that node. This is defined
        % once per N, gamma and is set global to save time.
        neighbours = zeros(n);
        neighbourCount = zeros(1, n);
        for ne = 1:N
            nbs = find(W(ne,:));
            neighbourCount(ne) = length(nbs);
            neighbours(ne,1:neighbourCount(ne)) = nbs;
        end
        rhoidx = 0;
        for rho = rhospan
            rhoidx = rhoidx + 1;
            iteridx = 0;
            for iter = iterations
                tstart = tic;
                iteridx = iteridx+1;
                % Randomly populate the nodes with mean rho
                pop = randomPop(rho, n);
                totalpop = sum(pop);
                % Set ~10% of the initial population to be infected and
                % ~90% to be susceptible.
                r = rand(1, n);
                I0 = pop.*(r>0.9);
                S0 = pop-I0;
                pops = [S0;I0];
                
                %I = zeros(tmax, n);
                % Calculate the relative population of I in relation to
                % rho.
                Irelative(:) = 0;
                for t = 1:tmax
                    pops = myCombinedSpread(dt, pops, W, status);
                    Irelative(t) = sum(pops(2,:))/totalpop;
                    if(Irelative(t) == 0)
                        break
                    end
%                     if(mod(t,101) == 0) && (abs(mean(Irelative((t-50):t))-mean(Irelative((t-100):(t-50)))) < tol)
%                         break
%                     end
                end
                span = (t-50):t;
                % Saving Irelative and total population.
                if(s)
                    infectiousFile = sprintf('INFECTEDrel-N%d-gamma%.1f-rho%.1f-iter%d-%s.mat',...
                        N,gamma,rho,iter,filename)
                    popFile = sprintf('initPop-N%d-gamma%.1f-rho%.1f-iter%d-%s.mat',...
                        N,gamma,rho,iter,filename)
                    save(infectiousFile, 'Irelative')
                    save(popFile, 'totalpop')
                end
                if(Irelative(t) ~= 0)
                    avgIrel(iteridx, rhoidx, nidx) = mean(Irelative(span));
                end
                % Plotting Irelative vs. t
                plot(1:t, Irelative(1:t))
                hold on
                tend = toc(tstart)
            end
        end
    end
end
if(s)
    avgIrelfile = sprintf('averageIrelative-%s',filename)
    save(avgIrelfile, 'avgIrel')
end