function solve_fishery(param,weights)
% builds networks
% solves networks

Model=param.model;
pr=param.pr; %probability of successful reservation
pd=param.pd; %probability of development
pir=param.pir; %probability of being successfully influenced to not develop or reserve by reserved neighbours
pid=param.pid; %probability of being successfully influenced to develop or not reserve by developed neigbours
cost=param.cost;  %management cost
G1=param.G1; %fishery network weights
G2=param.G2; %fishery network no weights
BR=param.BR; %bspecies distributions

%generate solution files

if weights==0
    filename='ArtisanFisheryNoWeights80';
    G=G2;
else
    filename='ArtisanFisheryWeights80';
    G=G1;
end

for ipir=1:length(pir) %loop through the pir
    for ipid=1:length(pid) %loop through the pid
        buildSPUDDNetConsPlanFishery(G,BR,pr,pd,pir(ipir),pid(ipid),cost,filename,Model); % build SPUDD file and save BR
        solveNetConsPlanFishery(G,pr,pd,pir(ipir),pid(ipid),filename,Model); %solve
    end
end

end
