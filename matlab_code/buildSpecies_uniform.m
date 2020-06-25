function BR=buildSpecies_uniform(n,n_species)

BR=zeros(n,n_species); %assumes that the number of species is divisble by number of nodes
d=fix(n_species/n);
for i=0:n-1
    BR(i+1,i*d+1:(i+1)*d)=1;
end

end
