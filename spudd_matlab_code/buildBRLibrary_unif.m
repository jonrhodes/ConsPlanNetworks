function buildBRLibrary_unif(param)
% generates and saves uniform species distribution

BR.nsites=param.n_sites; % number of sites
BR.nspecies=param.n_species; % number of species

% generate uniform distribution
BR.mat=build_uniform(BR.nsites,BR.nspecies);
BR.name=['UniformBRSites',num2str(BR.nsites),'NumSp',num2str(BR.nspecies)];
filename=BR.name;
save(['../problems/ConsMDP/','BR_library/',filename,'.mat'],'BR');

end
