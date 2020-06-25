function G=getG(name,n)
% name can be {'line','island','star','wheel','islandstar','islandline','starline'}

if strcmp(name,'line')
    G=build_line(n);
elseif strcmp(name,'island')
    G=build_island(n);
elseif strcmp(name,'star')
    G=build_star(n);
elseif strcmp(name,'wheel')
    G=build_wheel(n);
elseif strcmp(name,'islandstar')
    G=build_islandstar(n);
elseif strcmp(name,'islandline')
    G=build_islandline(n);
elseif strcmp(name,'starline')
    G=build_starline(n);
else
    'error'
end

end
