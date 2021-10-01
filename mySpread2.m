function popOut = mySpread2(W, pop)

% move = zeros(size(W,1));
popOut = zeros(1,length(pop));

for n = 1:size(W,1)
    % Get all neighbours for node n
    neighbours = find(W(n,:));
    % Get amount of neighbours nc
    nc = length(neighbours);
    % Generate a random vector with indices corresponding to neighbours
    % which will get an individual from node n
    indicies = neighbours(randi([1 nc],1,pop(n)));
    % Update popOut which gives a full description of how all
%     % individuals move.
    for i = 1:pop(n)
%         move(n,indicies(i)) = move(n,indicies(i)) + 1;
        popOut(indicies(i)) = popOut(indicies(i)) + 1;
    end
end

% Wout = sparse(move);