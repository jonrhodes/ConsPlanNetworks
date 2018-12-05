function G=build_line(n)

G=zeros(n);
for i=1:n-1
    G(i,i+1)=1;
    G(i+1,i)=1;
end

end


