function pop = randomPop(rho, n)

tol = 1e-5;
pop = round(0.5*rho*randi([-1 1], 1, n) + rho);

while(abs(sum(pop)-n*rho) > tol)
    if(sum(pop)-n*rho > tol)
        viable = find(pop > 0);
        l = length(viable);
        idx = viable(randi(l));
        pop(idx) = pop(idx)-1;
    end
    if(n*rho-sum(pop) > tol)
        idx = randi(n);
        pop(idx) = pop(idx) + 1;
    end
end

% while(sum(pop)-n*rho > 0)
%     viable = find(pop > 0);
%     l = length(viable);
%     idx = viable(randi(l));
%     pop(idx) = pop(idx)-1;
% end
% while(sum(pop)-n*rho < 0)
%     idx = randi(n);
%     pop(idx) = pop(idx) + 1;
% end

