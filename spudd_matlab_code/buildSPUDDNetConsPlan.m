%% 13 feb 2013
% iadine.chades@csiro.au
% paper: Conservation planning across time and actors

% This is the main file to call to define the problem (FMDP)

function buildSPUDDNetConsPlan(G,BR,pr,pd,pir,pid,cost,Shape,filename,Model)

% build factored MDP
ddFMDP=mdpConsNetworkSIR(G,pr,pd,pir,pid,cost,BR,Model);
% save factored MDP
printSIRMDPspuddFormat(['../problems/ConsMDP/',Shape,'/','Model',num2str(Model),filename,'pr',num2str(pr),'pd',num2str(pd),'pir',num2str(pir),'pid',num2str(pid),'.dat'], ddFMDP);

end
