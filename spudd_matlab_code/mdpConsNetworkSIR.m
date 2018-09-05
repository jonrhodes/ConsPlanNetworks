function ddMDP=mdpConsNetworkSIR(G,pr,pd,pir,pid,cost,BR,Model)
% Author: Iadine Chades
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
% G: adjacency matrix; pr, pd, pid, pir, BR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
% ddMDP: mdp struct (factored MDP) encoded with algebraic decision diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments:
%
% Generates a POMDP struct (factored POMDP) of the weeds network problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(G)>20
    error('The adjacency graph G must satisfy n < 20');
end

n=length(G);

parentParcels = cell(n,1); % list of neighbors and weight of the interactions
for i = 1:n
    count=1;
    for j=1:n
        
        if j~=i && G(i,j)~=0
            parentParcels{j} = [parentParcels{j}, i];
            count=count+1;
        end
        % empty indicates no parent
    end
end

if Model==1 %reserved and developed both influence development
    ddMDP=fillddMDPNetSIR1(pr,pd,pir,pid,cost,BR,parentParcels);
elseif Model==2 %reserved and developed both influence reservation
    ddMDP=fillddMDPNetSIR2(pr,pd,pir,pid,cost,BR,parentParcels);
elseif Model==3 %reserved influences reservation and developed influences development
    ddMDP=fillddMDPNetSIR3(pr,pd,pir,pid,cost,BR,parentParcels);
else
    error('Method must equal 1, 2, or 3')
end

end