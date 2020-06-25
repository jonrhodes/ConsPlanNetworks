%% recursive Decision Diagrams building function
% This function builds the CPT using a decision diagram (DD) library
% INPUT:
% nliste : list of nodes connected to the node Id (list of variables to be
% considered)
% Id : node or variable for which we are writing the CPT. The state of Id
% depends on the states of its list of connected nodes.
% pr : probability that we are computing for each branch/nodes of the tree
% col: probability of colonisation from connected nodes to Id
% flag=1 if we are looking for P(Id'=reserved|All other nodes)
%   = 0 if we are looking for P(Id'=available|All other nodes)
%   = 2 if we are looking for P(Id'=reserved|All other nodes)

function tree_DD = BuildDDMarginalRew(sitelist,BR,new_siteR,sitesR)

if isempty(sitelist)
    % Leaf!

    tree_DD = DDleaf.myNew(sum(max(sitesR,new_siteR)-sitesR)); % on fait la somme - et on renvoit la feuille

else
    % node
    %%% BUG HERE
    longueur = size(sitelist,2);
    Id=sitelist(1);

    tree_DD = DDnode.myNew(Id, ...
        [BuildDDMarginalRew(sitelist(2:1:longueur),BR,new_siteR,sitesR), ...    % Available sitesR is unchanged
        BuildDDMarginalRew(sitelist(2:1:longueur),BR,new_siteR,max(sitesR,BR(sitelist(1),:))), ...     % Reserved we might get some new species
        BuildDDMarginalRew(sitelist(2:1:longueur),BR,new_siteR,sitesR)]);       % Dev sitesR is unchanged
end
