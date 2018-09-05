function G=build_islandstar(n)

G=zeros(n);

%get size of islands and stars
n_island=round(n/2);
n_star=n-n_island;

%build island
G_island=build_island(n_island);

%build star
G_star=build_star(n_star);

%populate G
G(1:n_island,1:n_island)=G_island;
G((n_island+1):n,(n_island+1):n)=G_star;

%join star and island
G(n,n_island)=1;
G(n_island,n)=1;

end
