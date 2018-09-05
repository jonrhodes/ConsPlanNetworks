function [graphM,terminalAct]=build_policy_tree(filename)
% input:

% filename is a policy dot file from spudd
% output:
% graphM is a (m,4) matrix representing the policy graph with m the amount
% of nodes, and corresponding links, graphM(node,1:3) [avail,res,dev]
% graphM(node,4)= the # corresponding site or -1 if terminal node
% if terminal node, actions can be found in terminalAct

[terminalAct, sitesID]=find_node(filename);
[yg,Clink]=find_link(filename);

[n foo]=size(terminalAct);
[foo m]=size(Clink);

graphM=zeros(n,4);

% go through each link and populate the matrix
for i=1:m
    inode=str2num(char(Clink{i}(1)))+1;
    next_node=str2num(char(Clink{i}(2)))+1;
    str_link=char(Clink{i}(3));
    
    if strcmp(str_link,'avail')
        graphM(inode,1)= next_node;
    elseif strcmp(str_link,'res')
        graphM(inode,2)= next_node;
    elseif strcmp(str_link,'dev')
        graphM(inode,3)= next_node;
    else
        'error in build_policy_tree'
    end
end
% we have the links we need the actions
graphM(:,4)=sitesID(:,2);
end
