% recursive Decision Diagrams building function
% This function builds the CPT using a decision diagram (DD) library
% author: iadine chades
% INPUT:
% nliste : list of nodes connected to the node Id (list of variables to be
% considered)
% wliste : list of node weights connected to the node Id (list of variable weights to be
% considered)
% Id : node or variable for which we are writing the CPT. The state of Id
% depends on the states of its list of connected nodes.
% pr : probability that we are computing for each branch/nodes of the tree
% col: probability of colonisation from connected nodes to Id
% flag=1 if we are looking for P(Id'=reserved|All other nodes)
%   = 0 if we are looking for P(Id'=available|All other nodes)
%   = 2 if we are looking for P(Id'=reserved|All other nodes)

function tree_DD = BuildDDSIR3Fishery(nliste,wliste,Id,pr,pd,pir,pid,nr,nd,flag,neighbours)
precision_n=1000;
if isempty(nliste)
    if (neighbours == 0)
        % Leaf!
        pR=round(pr*precision_n)/precision_n;
        pD=abs(max(0,round((1-pR)*(pd*precision_n))/precision_n));
        pA=abs(max(0,round((1-pD-pR)*precision_n)/precision_n));        
    else
        % Leaf!
        pR=round((pr+(1-pr)*(1-prod(1-(wliste*pir.*nr))))*precision_n)/precision_n;
        pD=abs(max(0,round((1-pR)*(pd+(1-pd)*(1-prod(1-(wliste*pid.*nd))))*precision_n)/precision_n));
        pA=abs(max(0,round((1-pD-pR)*precision_n)/precision_n));
    end
        
    if flag == 1 % P(Id'= reserved| ...
        tree_DD = DDleaf.myNew(pR);
    elseif flag == 2 % P(Id'= developed| ...
        tree_DD = DDleaf.myNew(pD);
    elseif flag==0 % P(Id'= Avail| ...
        tree_DD =DDleaf.myNew(pA);
    end
else % neigbours
    % node
    longeur = size(nliste,2);
    if longeur == 1
        nrR = nr;
        nrR(size(wliste,2)) = 1;
        ndD = nd;
        ndD(size(wliste,2)) = 1;
    else
        nrR = nr;
        nrR(size(wliste,2) - size(nliste(2:1:longeur),2)) = 1;
        ndD = nd;
        ndD(size(wliste,2) - size(nliste(2:1:longeur),2)) = 1;
    end
    tree_DD = DDnode.myNew(nliste(1), ...
        [BuildDDSIR3Fishery(nliste(2:1:longeur),wliste,Id,pr,pd,pir,pid,nr,nd,flag,neighbours), ...    % Available
        BuildDDSIR3Fishery(nliste(2:1:longeur),wliste,Id,pr,pd,pir,pid,nrR,nd,flag,neighbours), ...     % Reserved
        BuildDDSIR3Fishery(nliste(2:1:longeur),wliste,Id,pr,pd,pir,pid,nr,ndD,flag,neighbours)]);       % Dev
end
