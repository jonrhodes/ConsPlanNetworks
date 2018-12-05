function G=build_island(n)

G=build_line(n);
G(1,end)=1;
G(end,1)=1;

end

