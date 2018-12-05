function G=build_wheel(n)

G=zeros(n);
% the last node is the central node;
G(n,1:n-1)=1;
G(1:n-1,n)=1;

Gnew=build_line(n-1);
Gnew(1,end)=1;
Gnew(end,1)=1;

G(1:n-1,1:n-1)=Gnew;

end
