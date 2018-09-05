function buildBRLibrary(param)
% generates and saves species distribution for a given level of nestedness

BR.nsites=param.n_sites; % number of sites
BR.nspecies=param.n_species; % number of species
BR.pocc=param.p_occurence; % probability of occurence at a site
BR.nested=param.nested; %nestedness
n_rep = param.n_rep; % number of replicates of the BR matrix

for i=1:n_rep
    BR.rep=i;
    BR.mat=buildSpecies(BR.nsites,BR.nspecies,BR.pocc,BR.nested);
    BR.name=['BRSites',num2str(BR.nsites),'NumSp',num2str(BR.nspecies),'Occ',num2str(BR.pocc),'Nested',num2str(BR.nested),'Rep',num2str(i)];
    filename=BR.name;
    save(['../problems/ConsMDP/','BR_library/',filename,'.mat'],'BR');
end

end
