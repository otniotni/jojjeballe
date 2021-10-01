% Function networkgenerator

function net = networkgenerator(V, gamma, savefile)


xmin = 2;
xmax = sqrt(V);
net = zeros(V);
k = zeros(1, V);

powerdist = @(r) ((xmax.^(gamma+1)-xmin.^(gamma+1)).*r+xmin.^(gamma+1)).^(1/(1+gamma));

index = 1;
for node = 1:V
    k(node) = round(powerdist(rand));
    kc(index:index+k(node)) = node;
    index = index + k(node) + 1;
end

while ~isempty(kc)
    indices = length(kc);
    twoNodes = randi([1 indices], 1, 2);
    net(kc(twoNodes(1)),kc(twoNodes(2))) = 1;
    net(kc(twoNodes(2)),kc(twoNodes(1))) = 1;
    kc(twoNodes) = [];
end
net = sparse(net);
if ~isempty(savefile)
    save(sprintf('C:\\Users\\Otto\\Documents\\MODSIM\\Lab 2\\powernetsByOtto\\%s',savefile), 'net')
end




