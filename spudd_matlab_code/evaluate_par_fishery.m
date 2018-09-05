function evaluate_par_fishery(param)
% evaluate performance of each method under different scenario
% data from simulations are saved in Results/

Model=param.model;
pr=param.pr; %probability of successful reservation
pd=param.pd; %probability of development
pir=param.pir; %probability of being successfully influenced to not develop or reserve by reserved neighbours
pid=param.pid; %probability of being successfully influenced to develop or not reserve by developed neigbours
cost=param.cost;  %management cost
G1=param.G1; %fishery network weights
G2=param.G2; %fishery network no weights
BR=param.BR; %bspecies distributions
n_simu=param.n_simu; %number of simulations for each network, BR and solution methods

ExpR=zeros(length(pir),length(pid),1,2); %the results _ expected number of species saved
SDR=zeros(length(pir),length(pid),1,2);  %the results _ sd of expected number of species saved
ExpRc=zeros(length(pir),length(pid),1,2); %the results _ expected number of species saved including cost
SDRc=zeros(length(pir),length(pid),1,2);  %the results _ sd of expected number of species saved including cost

%SpecPerc = '95';
%SpecPerc = '80';
SpecPerc = '50';

parpool(30);

%evaluate for G1 - with differences in tie strength

for ipir=1:length(pir) %loop over the pir
    for ipid=1:length(pid) %loop over pid
        for meth=1:2 %loop over the methods - in this case it is just CPSN and SN
            %which methods?
            if meth==1 %conservation planning and social network
                filename=['ArtisanFisheryWeights',SpecPerc,'Model3pr0.2pd0.2','pir',num2str(pir(ipir)),'pid',num2str(pid(ipid))];
                pol_file=['../ConsMDP/fishery/',filename,'-OPTpolicy.dot'];
                [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                [ExpR(ipir,ipid,1,meth),SDR(ipir,ipid,1,meth),ExpRc(ipir,ipid,1,meth),SDRc(ipir,ipid,1,meth)]=run_policy_fishery(G1,BR,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
            elseif meth==2 %conservation planning - apply 'pid=0 and pir=0' optimal strategy
                filename=['ArtisanFisheryWeights',SpecPerc,'Model3pr0.2pd0.2','pir',num2str(0),'pid',num2str(0)];
                pol_file=['../ConsMDP/fishery/',filename,'-OPTpolicy.dot'];
                [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                [ExpR(ipir,ipid,1,meth),SDR(ipir,ipid,1,meth),ExpRc(ipir,ipid,1,meth),SDRc(ipir,ipid,1,meth)]=run_policy_fishery(G1,BR,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
            else
                error('Solution type does not exist');
            end
        end
    end
end


save(['../ConsMDP/','results/ExpR_FisheryWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpR');
save(['../ConsMDP/','results/SDR_FisheryWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDR');
save(['../ConsMDP/','results/ExpRc_FisheryWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpRc');
save(['../ConsMDP/','results/SDRc_FisheryWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDRc');

%evaluate for G2 - with no differences in tie strength

for ipir=1:length(pir) %loop over the pir
    for ipid=1:length(pid) %loop over pid
        for meth=1:2 %loop over the methods - in this case it is just CPSN and SN
            %which methods?
            if meth==1 %conservation planning and social network
                filename=['ArtisanFisheryNoWeights',SpecPerc,'Model3pr0.2pd0.2','pir',num2str(pir(ipir)),'pid',num2str(pid(ipid))];
                pol_file=['../ConsMDP/fishery/',filename,'-OPTpolicy.dot'];
                [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                [ExpR(ipir,ipid,1,meth),SDR(ipir,ipid,1,meth),ExpRc(ipir,ipid,1,meth),SDRc(ipir,ipid,1,meth)]=run_policy_fishery(G2,BR,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
            elseif meth==2 %conservation planning - apply 'pid=0 and pir=0' optimal strategy
                filename=['ArtisanFisheryNoWeights',SpecPerc,'Model3pr0.2pd0.2','pir',num2str(0),'pid',num2str(0)];
                pol_file=['../ConsMDP/fishery/',filename,'-OPTpolicy.dot'];
                [MG,MA]=build_policy_tree(pol_file); %MG and MA are the optimal policy tree;
                [ExpR(ipir,ipid,1,meth),SDR(ipir,ipid,1,meth),ExpRc(ipir,ipid,1,meth),SDRc(ipir,ipid,1,meth)]=run_policy_fishery(G2,BR,MG,MA,pr,pd,pir(ipir),pid(ipid),n_simu,Model); % call function defined at the bottom
            else
                error('Solution type does not exist');
            end
        end
    end
end

save(['../ConsMDP/','results/ExpR_FisheryNoWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpR');
save(['../ConsMDP/','results/SDR_FisheryNoWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDR');
save(['../ConsMDP/','results/ExpRc_FisheryNoWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'ExpRc');
save(['../ConsMDP/','results/SDRc_FisheryNoWeights',SpecPerc,'Model',num2str(3),'pr',num2str(pr),'pd',num2str(pd),'.mat'],'SDRc');

delete(gcp);

end

function [Aver,SD,Averc,SDc] =run_policy_fishery(G,BR,MG,MA,pr,pd,pir,pid,n_simu,model) % call function defined at the bottom

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
            %Snew2 = simNetwork1(Snew,G,Act,pr,pd,pir,pid);
            error('Model type not supported');
        elseif model==2
            %Snew2 = simNetwork2(Snew,G,Act,pr,pd,pir,pid);
            error('Model type not supported');
        elseif model==3
            Snew2 = simNetwork3_fishery(Snew,G,Act,pr,pd,pir,pid);
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
end

Aver=mean(R);
SD=std(R)/sqrt(n_simu);
Averc=mean(C);
SDc=std(C)/sqrt(n_simu);

end
