% Function popOut distributes the population in pop according to the links
% in W. This function treats the population in pop as two separate
% species: infected and susceptible individuals.

function popOut = mySpread2mod(W, pop)
global neighbours neighbourCount status
species = size(pop,1);
popOut = zeros(size(pop));
nspan = 1:size(W,1);
if(status == 5) 
    Irel = pop(2,:)./sum(pop);
    threshold = 0.7;
%     noQuarantine = find(Irel<0.5);
    % The following actions puts all nodes with an Irelative over
    % threshold in quarantine, meaning that no one can leave this node.
    % Individuals can however enter the node.
    nspan = nspan(Irel<threshold);
    popOut(:,Irel >= threshold) = pop(:,Irel >= threshold); 
end

if(species == 2)
    for n = nspan
        nc = neighbourCount(n);
        neighs = neighbours(n,1:nc);
        % Generate a random vector with indices corresponding to neighbours
        % which will get an individual from node n
        indicies1 = neighs(randi([1 nc],1,pop(1,n)));
        indicies2 = neighs(randi([1 nc],1,pop(2,n)));
        % Update popOut which gives a full description of how all
        for i = 1:max(length(indicies1),length(indicies2))
            if(i <= length(indicies1))
                popOut(1,indicies1(i)) = popOut(1,indicies1(i)) + 1;
            end
            if(i <= length(indicies2))
                popOut(2,indicies2(i)) = popOut(2,indicies2(i)) + 1;
            end
        end
    end
elseif(species == 1)
    for n = 1:size(W,1)
        nc = neighbourCount(n);
        neighs = neighbours(n,1:nc);
        % Generate a random vector with indices corresponding to neighbours
        % which will get an individual from node n
        indicies1 = neighs(randi([1 nc],1,pop(1,n)));
        %indicies2 = neighs(randi([1 nc],1,pop(2,n)));
        % Update popOut which gives a full description of how all
        for i = 1:length(indicies1)
            popOut(1,indicies1(i)) = popOut(1,indicies1(i)) + 1;
        end
    end
end
