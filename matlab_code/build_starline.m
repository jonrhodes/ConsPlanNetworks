function G=build_starline(n)

G=zeros(n);

%get size of stars and lines
n_star=round(n/2);
n_line=n-n_star;

%build star
G_star=build_star(n_star);

%build line
G_line=build_line(n_line);

%populate G
G(1:n_star,1:n_star)=G_star;
G((n_star+1):n,(n_star+1):n)=G_line;

%join star and line
G(n,n_star)=1;
G(n_star,n)=1;

end
