function main_fishery(Model)
%solves fishery case study
%Model=1 (reserved and developed both influence development) or 2
%(reserved and developed both influence reservation) or 3 (reserved
%influences reservation and developed influences development)
%nSites=number of nodes

tic

%get social network 
G=build_fishery;
%standardise tie strength based on the maximum tie strength
G1=G/max(max(G));
%ignore weights but keep the same average tie strength
G2=sign(G1)*(sum(sum(G1))/sum(sum(sign(G1))));

%build the fishery species matrix
%BR=buildFisherySpecies('fishery_species95.csv');
%BR=buildFisherySpecies('fishery_species80.csv');
BR=buildFisherySpecies('fishery_species50.csv');


%define parameter values
param.G1=G1;
param.G2=G2;
param.BR=BR;
param.model=Model;
param.pr=0.2; %probability of successful reservation
param.pd=0.2; %probability of development
param.pir=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %probability of being successfully influenced to not develop or reserve by reserved neighbours
param.pid=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %probability of being successfully influenced to develop or not reserve by developed neigbours
param.cost=-0.5; %management cost to incorporate
param.n_simu = 100000; %number of simulations for each network, BR and solution methods

%solve with and without weights
solve_fishery(param,0);
solve_fishery(param,1);

toc

end
