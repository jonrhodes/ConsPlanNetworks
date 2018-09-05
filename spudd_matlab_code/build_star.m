
function G=build_star(n)

G=zeros(n);
% the last node is the central node;
G(n,1:n-1)=1;
G(1:n-1,n)=1;

end
