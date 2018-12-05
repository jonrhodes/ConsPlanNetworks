%% recursive Decision Diagrams building function
% This function builds the CPT using a decision diagram (DD) library
% author: iadine chades
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

function tree_DD = BuildDDSIR1(nliste,Id,pr,pd,pir,pid,nr,nd,flag)
precision_n=1000;
if isempty(nliste)
    % Leaf!
    pD=round((1-pr)*(pd+(1-pd)*(1-(1-pid)^nd))*((1-pir)^nr)*precision_n)/precision_n;
    pR=pr;
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
        [BuildDDSIR1(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr,nd,flag), ...    % Available
        BuildDDSIR1(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr+1,nd,flag), ...     % Reserved
        BuildDDSIR1(nliste(2:1:longueur),Id,pr,pd,pir,pid,nr,nd+1,flag)]);       % Dev
end
