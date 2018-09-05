function simNetConsPlan()
%% main file to call
% builds networks
% solves networks
%%

%ShapeC={'line','island','star','cluster'};
ShapeC={'line'};        % cluster is too long to run for now
%BiodC={'uniform','Rand'};
BiodC={'Rand'};
% probabilities
pr= 0.3;    % prob of successful reservation
pd= 0.1;    % prob of development
%pid= [0., 0.1, 0.2, 0.3, 0.4, 0.5 ];   % prob of being successfully influenced to develop a site by developed neighbours
%pidname={'pzero', 'pone','ptwo','pthree','pfour','pfive'};   % has to match the pid vector and directory name
pid= [0.1,0.3 0.5];   % prob of being successfully influenced to develop a site by developed neighbours
pidname={'pone','pthree','pfive'};   % has to match the pid vector and directory name


%pid= [0.5 ];   % prob of being successfully influenced to develop a site by developed neighbours
%pidname={'pfive'};   % has to match the pid vector and directory name


pir= pid;   % prob of being successfully influenced to not develop a site by reserved neigbours
cost=-0.5;  % management cost

% number of nodes
n=10;
maxr=80;

% generate all solution files for combinaison of shape and BR
for i=1:length(ShapeC)
    % build the net - can also call getG(..)
    Sshape=ShapeC{i};
    if strcmp(Sshape,'line')
        G=build_line(n);
    elseif strcmp(Sshape,'star')
        G=build_star(n);
    elseif strcmp(Sshape,'cluster')
        G=build_cluster(n);
    elseif strcmp(Sshape,'island')
        G=build_island(n);
    else
        Sshape='rand';  % unused atm
        G=rand(n);
    end
    
    replicate_id=1;
    % buid the BR
    for j=1:length(BiodC)   % loop through the BR
        Biod=BiodC{j}
        if strcmp(Biod,'uniform')
            BR=load(['../problems/ConsMDP/BR_library/uniformBR.mat'],'BR');  % previously saved
            for ipid=1:length(pid) %loop through the pid
                buildSPUDDNetConsPlan(G,BR.BR,pr,pd,pid(ipid),pidname{ipid},pir(ipid),cost,Sshape,Biod,replicate_id); % build SPUDD file and save BR
                solveNetConsPlan(G,BR.BR,Sshape,Biod,replicate_id,pidname{ipid});    % solve
            end
        else
            Biod='rand';
            for rep_id=1:10:maxr
                BR=load(['../problems/ConsMDP/BR_library/BR',num2str(rep_id),'.mat']);  % previously built using buildBRLibrary.m
                for ipid=1:length(pid)
                    buildSPUDDNetConsPlan(G,BR.x,pr,pd,pid(ipid),pidname{ipid},pir(ipid),cost,Sshape,Biod,rep_id);
                    solveNetConsPlan(G,BR.x,Sshape,Biod,rep_id,pidname{ipid});
                end
            end
        end
    end
    
    
end
end



