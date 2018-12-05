function solve_all(param)
% builds networks
% solves networks

Model=param.model;
n_species=param.n_species; %number of species
Nested=param.nested; %nestedness (0 = perfectly nested to 0.5 = random)
n_sites=param.n_sites; %number of sites
Shape=param.shape; %shape of network
pr=param.pr; %probability of successful reservation
pd=param.pd; %probability of development
pir=param.pir; %probability of being successfully influenced to not develop or reserve by reserved neighbours
pid=param.pid; %probability of being successfully influenced to develop or not reserve by developed neigbours
cost=param.cost;  %management cost
n_rep=param.n_rep; %number of bBR replicates to solve
p_occurrence=param.p_occurrence; %probability of occupancy

%generate all solution files for each BR

%build the network given the shape and the number of sites
G=getG(Shape,n_sites);
%build the BR
for rep_id=1:n_rep %loop through the replicates to be solved
    filename=['BRSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'Rep',num2str(rep_id)];
    load(['../problems/ConsMDP/','BR_library/',filename,'.mat']); %previously built using buildBRLibrary.m
    for ipir=1:length(pir) %loop through the pir
        for ipid=1:length(pid) %loop through the pid
            ['Solve Other: ','Pir = ',num2str(pir(ipir)),' Pid = ',num2str(pid(ipid))]
            buildSPUDDNetConsPlan(G,BR.mat,pr,pd,pir(ipir),pid(ipid),cost,Shape,filename,Model); % build SPUDD file and save BR
            solveNetConsPlan(G,Shape,pr,pd,pir(ipir),pid(ipid),filename,Model); %solve
        end
    end
end

end
