%% recursive Decision Diagrams building function

function tree_DD = BuildDDRew(sitelist,BR,sitesR)

if isempty(sitelist)
    % Leaf!
    
    tree_DD = DDleaf.myNew(sum(sitesR)); % on fait la somme - et on renvoit la feuille
    
else
    longueur = size(sitelist,2);
    Id=sitelist(1);
    
    tree_DD = DDnode.myNew(Id, ...
        [BuildDDRew(sitelist(2:1:longueur),BR,sitesR), ...    % Available sitesR is unchanged
        BuildDDRew(sitelist(2:1:longueur),BR,max(sitesR,BR(sitelist(1),:))), ...     % Reserved we might get some new species
        BuildDDRew(sitelist(2:1:longueur),BR,sitesR)]);       % Dev sitesR is unchanged
end
