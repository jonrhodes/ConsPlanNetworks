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

function tree_DD = BuildDDSIR2(nliste,Id,pr,pd,pir,pid,nr,nd,flag)
precision_n=1000;
if isempty(nliste)
    % Leaf!
    pR=round((pr+(1-pr)*(1-(1-pir)^nr))*((1-pid)^nd)*precision_n)/precision_n;
    pD=abs(max(0,round((1-pR)*pd*precision_n)/precision_n));
    pA=abs(max(0,round((1-pD-pR)*precision_n)/precision_n));

    if flag == 1 % P(Id'= reserved| ...
        tree_DD = DDleaf.myNew(pR);
    elseif flag == 2 % P(Id'= developed| ...
        tree_DD = DDleaf.myNew(pD);
    elseif flag==0 % P(Id'= Avail| ...
        tree_DD =DDleaf.myNew(pA);
    end
else % neigbours
    % node
    longueur = size(nliste,2);
    tree_DD = DDnode.myNew(nliste(1), ...
        [BuildDDSIR2(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr,nd,flag), ...    % Available
        BuildDDSIR2(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr+1,nd,flag), ...     % Reserved
        BuildDDSIR2(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr,nd+1,flag)]);       % Dev
end
