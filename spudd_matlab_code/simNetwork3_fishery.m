function [SNew] = simNetwork3_fishery(S,N,Act,pr,pd,pir,pid)
%outputs a new state SNew given and initial state (S)
%an action (A) and the parameters pr, pd, pir, pid
%S = n x 1 matrix describing the state with 0 = available, 1 = reserved,
%and 2 = developed
%N = n x n matrix describing the network
%Act = number between 1 and n indicating which site is reserved (the action)
%pr = probability of reservation given site reserved
%pd = probability of development given site available
%pir = probability site infected from connection to reserved site
%pid = probability site is infected from connection to developed site

precision_n=1000;

%find available sites not being reserved
Avail=find(S==0);
Avail=Avail(find(Avail~=Act));

%find reserved sites
Res=find(S==1);

%find developed sites
Dev=find(S==2);

%simulate transitions to available/reserved/developed for available sites
%being reserved

%simulate transitions to available/developed for available sites not being reserved

if ~isempty(Avail)    
    %find the connections to reserved and developed sites
    Nr=N(Avail,Res);
    Nd=N(Avail,Dev);
    
    %get transtion probabilies when not reserving
    pD=round((pd+(1-pd)*(1-prod(1-(pid.*Nd),2)))*precision_n)/precision_n;
       
    %simulate
    AvailtoAD=binornd(1,pD);
    S(Avail)=S(Avail)+(2.*AvailtoAD);    
end
    
%find number of connections to reserved and developed sites
if Act~=0  % Act==0 is the do nothing action
    if S(Act)~=0
        'action aborted'
    else
        Nr=N(Act,Res);
        Nd=N(Act,Dev);

        %get transtion probabilities
        pR=round((pr+(1-pr)*(1-prod(1-(pir.*Nr),2)))*precision_n)/precision_n;
        pD=abs(max(0,round((1-pR).*(pd+(1-pd)*(1-prod(1-(pid.*Nd),2)))*precision_n)/precision_n));
        pA=abs(max(0,round((1-pD-pR)*precision_n)/precision_n));

        %simulate
        AvailtoDRA=mnrnd(1,[pD pR pA]);
        %store transitions for available site being reserved
        S(Act)=S(Act)+(2.*AvailtoDRA(1))+AvailtoDRA(2);
    end
end

SNew=S;

end