function G=build_islandline(n)

G=zeros(n);

%get size of islands and stars
n_island=round(n/2);
n_line=n-n_island;

%build island
G_island=build_island(n_island);

%build line
G_line=build_line(n_line);

%populate G
G(1:n_island,1:n_island)=G_island;
G((n_island+1):n,(n_island+1):n)=G_line;

%join island and line
G(n,n_island)=1;
G(n_island,n)=1;

end
