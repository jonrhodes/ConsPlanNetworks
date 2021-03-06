% This is the main file to call to define and solve the problem (FMDP)

function solveNetConsPlanFishery(G,pr,pd,pir,pid,filename,Model)

Consdir=['../problems/ConsMDP/fishery/'];

% call SPUDD solver

filename=[filename,'Model',num2str(Model),'pr',num2str(pr),'pd',num2str(pd),'pir',num2str(pir),'pid',num2str(pid)];

cmd=['$HOME/spudd-3.6.2/Spudd/bin/linux/Spudd ',Consdir,filename,'.dat -dd -o ',Consdir,filename];
system(cmd);

% generate policy solution PDF
cmd=['dot -Tpdf ',Consdir,filename,'-OPTpolicy.dot -o ',Consdir, filename,'_Sol.pdf'];
system(cmd);

% plot studied graph
graph_to_dot(G,'filename',[Consdir,filename,'_Net.dot'],'directed',0);
cmd=['dot -Tpdf ',Consdir,filename,'_Net.dot -o ',Consdir, filename,'_Net.pdf'];
system(cmd);

end
