function evaluate_par(param)
% evaluate performance of each method under different scenario
% data from simulations are saved in Results/

Model=param.model;
n_species=param.n_species; %number of species
Nested=param.nested;
n_sites=param.n_sites; %number of sites
Shape=param.shape; %shape of network
pr=param.pr; %probability of successful reservation
pd=param.pd; %probability of development
pir=param.pir; %probability of being successfully influenced to not develop or reserve by reserved neighbours
pid=param.pid; %probability of being successfully influenced to develop or not reserve by developed neigbours
cost=param.cost;  %management cost
n_rep=param.n_rep; %number of BR replicates to evaluate
p_occurrence=param.p_occurrence; %probability of occupancy
n_simu=param.n_simu; %number of simulations for each network, BR and solution methods
sol_methods=param.sol_methods; %CPSN= us; CP= no social network (pid=pir=0); SN= only SN = uniform biodiversity distribution

ExpR=zeros(length(pir),length(pid),n_rep,length(sol_methods)); %the results _ expected number of species saved
SDR=zeros(length(pir),length(pid),n_rep,length(sol_methods));  %the results _ sd of expected number of species saved
ExpRc=zeros(length(pir),length(pid),n_rep,length(sol_methods)); %the results _ expected number of species saved including cost
SDRc=zeros(length(pir),length(pid),n_rep,length(sol_methods));  %the results _ sd of expected number of species saved including cost

G=getG(Shape,n_sites);

parpool(30);

for rep_id=1:n_rep %loop over the replicates
    
    % get the corresponding BR
    filename=['BRSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'Rep',num2str(rep_id)];
    load(['../ConsMDP/','BR_library/',filename,'.mat']); %previously built using buildBRLibrary.m
    
    for ipir=1:length(pir) %loop over the pir
        for ipid=1:length(pid) %loop over pid
            for meth=1:length(sol_methods) %loop over the methods
                %['Evaluate: ','Pir = ',num2str(pir(ipir)),' Pid = ',num2str(pid(ipid)),' Method = ',sol_methods{meth}]                
                %which methods?
                if strcmp(sol_methods{meth},'CPSN') %conservation planning and social network
                    filename=['BRSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'Rep',num2str(rep_id)];
                    filename=['Model',num2str(Model),filename,'pr',num2str(pr),'pd',num2str(pd),'pir',num2str(pir(ipir)),'pid',num2str(pid(ipid))];
                    pol_file=['../ConsMDP/',Shape,'/',filename,'-OPTpolicy.dot'];
                    [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                    [ExpR(ipir,ipid,rep_id,meth),SDR(ipir,ipid,rep_id,meth),ExpRc(ipir,ipid,rep_id,meth),SDRc(ipir,ipid,rep_id,meth)]=run_policy(G,BR.mat,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
                    
                elseif strcmp(sol_methods{meth},'CP') %conservation planning - apply 'pid=0 and pir=0' optimal strategy
                    filename=['BRSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'Rep',num2str(rep_id)];
                    filename=['Model',num2str(Model),filename,'pr',num2str(pr),'pd',num2str(pd),'pir',num2str(0),'pid',num2str(0)];
                    pol_file=['../ConsMDP/',Shape,'/',filename,'-OPTpolicy.dot'];
                    [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                    [ExpR(ipir,ipid,rep_id,meth),SDR(ipir,ipid,rep_id,meth),ExpRc(ipir,ipid,rep_id,meth),SDRc(ipir,ipid,rep_id,meth)]=run_policy(G,BR.mat,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
                    
                elseif strcmp(sol_methods{meth},'SN') %social network - the BR doesn't matter apply the 'uniformBR' optimal strategy
                    filename=['UniformBRSites',num2str(n_sites),'NumSp',num2str(n_species)];
                    filename=['Model',num2str(Model),filename,'pr',num2str(pr),'pd',num2str(pd),'pir',num2str(pir(ipir)),'pid',num2str(pid(ipid))];
                    pol_file=['../ConsMDP/',Shape,'/',filename,'-OPTpolicy.dot'];
                    [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                    [ExpR(ipir,ipid,rep_id,meth),SDR(ipir,ipid,rep_id,meth),ExpRc(ipir,ipid,rep_id,meth),SDRc(ipir,ipid,rep_id,meth)]=run_policy(G,BR.mat,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
                    
                elseif strcmp(sol_methods{meth},'RND') %random selection of actions
                    [ExpR(ipir,ipid,rep_id,meth),SDR(ipir,ipid,rep_id,meth),ExpRc(ipir,ipid,rep_id,meth),SDRc(ipir,ipid,rep_id,meth)]=run_policy_rnd(G,BR.mat,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
                               
                else
                    error('Solution type does not exist');
                end
            end
        end
    end
end

save(['../ConsMDP/','results/ExpR_',Shape,'Model',num2str(Model),'nSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpR');
save(['../ConsMDP/','results/SDR_',Shape,'Model',num2str(Model),'nSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDR');
save(['../ConsMDP/','results/ExpRc_',Shape,'Model',num2str(Model),'nSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpRc');
save(['../ConsMDP/','results/SDRc_',Shape,'Model',num2str(Model),'nSites',num2str(n_sites),'NumSp',num2str(n_species),'Occ',num2str(p_occurrence),'Nested',num2str(Nested),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDRc');

delete(gcp);

end

function [Aver,SD,Averc,SDc] =run_policy(G,BR,MG,MA,pr,pd,pir,pid,n_simu,model) % call function defined at the bottom

R=zeros(n_simu,1);
C=zeros(n_simu,1);
T=150;   %finite time horizon
n=length(G);
Disc = 0.96; %discount factor

parfor k=1:n_simu
    Snew=zeros(n,1); %all available
    RewardR=0;
    RewardC=0;
    for t=1:T
        Act=getAction(Snew+1,MG,MA);
        if model==1
            Snew2 = simNetwork1(Snew,G,Act,pr,pd,pir,pid);
        elseif model==2
            Snew2 = simNetwork2(Snew,G,Act,pr,pd,pir,pid);
        elseif model==3
            Snew2 = simNetwork3(Snew,G,Act,pr,pd,pir,pid);
        else
            error('Model type does not exist');
        end
                              
        if Act~=0
            Cost=0.5;
        else
            Cost=0;
        end
        
        Snew=Snew2;
        
        %get reward
        RewardR=RewardR+((get_reward(Snew,BR))*(Disc^t)); 
        RewardC=RewardC+((get_reward(Snew,BR)-Cost)*(Disc^t)); 
               
        if Snew~=0  %all sites are reserved or developed.
            break
        end
    end
    
    %add reward over infinite time horizon
    RewardR=RewardR+(get_reward(Snew,BR)/((1/Disc)-1))-get_reward(Snew,BR)*((1-(1/Disc)^t)/((1/Disc)-1));
    RewardC=RewardC+(get_reward(Snew,BR)/((1/Disc)-1))-get_reward(Snew,BR)*((1-(1/Disc)^t)/((1/Disc)-1));
        
    R(k)=RewardR;    
    C(k)=RewardC; 
    %[r,c]=size(BR);    
end
Aver=mean(R);
SD=std(R)/sqrt(n_simu);
Averc=mean(C);
SDc=std(C)/sqrt(n_simu);

end

function [Aver,SD,Averc,SDc] =run_policy_rnd(G,BR,pr,pd,pir,pid,n_simu,model) % call function defined at the bottom

R=zeros(n_simu,1);
C=zeros(n_simu,1);
T=150;   %finite time horizon
n=length(G);
Disc = 0.96; %discount factor

parfor k=1:n_simu
    Snew=zeros(n,1); %all available
    RewardR=0;
    RewardC=0;    
    for t=1:T
        Avail=find(Snew==0);
        Act=Avail(randsample(length(Avail),1));
        if model==1
            Snew2 = simNetwork1(Snew,G,Act,pr,pd,pir,pid);
        elseif model==2
            Snew2 = simNetwork2(Snew,G,Act,pr,pd,pir,pid);
        elseif model==3
            Snew2 = simNetwork3(Snew,G,Act,pr,pd,pir,pid);
        else
            error('Model type does not exist');
        end
        
        if Act~=0
            Cost=0.5;
        else
            Cost=0;
        end
        
        Snew=Snew2;
        
        %get reward
        RewardR=RewardR+((get_reward(Snew,BR))*(Disc^t)); 
        RewardC=RewardC+((get_reward(Snew,BR)-Cost)*(Disc^t)); 
               
        if Snew~=0  %all sites are reserved or developed.
            break
        end
    end
    
    %add reward over infinite time horizon
    RewardR=RewardR+(get_reward(Snew,BR)/((1/Disc)-1))-get_reward(Snew,BR)*((1-(1/Disc)^t)/((1/Disc)-1));
    RewardC=RewardC+(get_reward(Snew,BR)/((1/Disc)-1))-get_reward(Snew,BR)*((1-(1/Disc)^t)/((1/Disc)-1));
        
    R(k)=RewardR;
    C(k)=RewardC;
    %[r,c]=size(BR);    
end

Aver=mean(R);
SD=std(R)/sqrt(n_simu);
Averc=mean(C);
SDc=std(C)/sqrt(n_simu);

end
