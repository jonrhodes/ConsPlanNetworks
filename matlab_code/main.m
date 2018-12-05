function main(Shape,Model,nSites)
$solves and evaluates simulations
%Shape='line','island','star','wheel','islandstar','islandline','starline'
%note here 'island' is referred to as 'ring' in the associated manuscript
%Model = 1 (reserved and developed both influence development) or = 2
%(reserved and developed both influence reservation) or = 3 (reserved
%influences reservation and developed influences development)
%nSites=number of nodes

tic

%define parameter values
ns=[20]; %total number of species - note that number of species must be divisible by number of nodes
nst=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %the level of nestedness (0 = maximum nestedness [but can't = 0] to 0.5 = random to 1 = maximum endemism)
poc=[0.25]; %probability of occupancy for each species at each node (this is a measure of rarity or abundance);
param.n_sites=nSites; %number of nodes
param.shape=Shape;
param.model=Model;
param.pr=0.2; %probability of successful reservation
param.pd=0.2; %probability of development
param.pir=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %probability of being successfully influenced to not develop or reserve by reserved neighbours
param.pid=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %probability of being successfully influenced to develop or not reserve by developed neigbours
param.cost=-0.5; %management cost to incorporate
param.n_rep = 5; % replication of the BR matrix
param.n_simu = 1000; %number of simulations for each network, BR and solution methods
param.sol_methods={'CPSN','CP','SN','RND'}; % CPandSN=species and social network; CP=species only (pid=pir=0); SN=only social network (uniform species distribution); RND=random actions

% %build BR library
for i=1:length(ns)
     param.n_species=ns(i);
     %build library for uniform species distribution
     buildBRLibrary_unif(param);
     for j=1:length(nst)
         param.nested=nst(j);
         for k=1:length(poc)
             param.p_occurence=poc(k);
             %build library of stochastic species distributions of given
             %nestdness
             buildBRLibrary(param);
         end
     end
 end

%solve
for i=1:length(ns)
    param.n_species=ns(i);
    %solve
    solve_all_unif(param);
    for j=1:length(nst)
        param.nested=nst(j);
        for k=1:length(poc)
            param.p_occurrence=poc(k);
            %solve
            solve_all(param);
        end
    end
end

%evaluate
for i=1:length(ns)
    param.n_species=ns(i);
    for j=1:length(nst)
        param.nested=nst(j);
        for k=1:length(poc)
            param.p_occurrence=poc(k);
            %evaluate
            evaluate_par(param);
        end
    end
end

toc

end
