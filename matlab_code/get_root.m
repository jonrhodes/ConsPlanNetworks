function [site_r,node_r]=get_root(Mgraph)
% the root site only appears once on the 4th column of Mgraph
% need to check that.
% input: Mgraph the policy graph
% output: 
% site_r 
% node_r
[n,x]= size(Mgraph);

i=1;
okay=0;

while okay==0 && i<=n
    site_r=Mgraph(i,4);
    node_r=1;
    okay=1;
    for k=1:i-1
        if site_r==Mgraph(k,4)  % site_r is not the root
            okay=0;
            i=i+1;
            break  
        end  
    end
    for k=i+1:n
        if site_r==Mgraph(k,4)  % site_r is not the root
            okay=0;
            i=i+1;
            break  
        end  
    end
end

end
