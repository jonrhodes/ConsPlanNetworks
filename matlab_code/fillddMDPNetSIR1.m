function ddMDP=fillddMDPNetSIR1(pr,pd,pir,pid,cost,BR,parentParcels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Iadine Chades  xx SIR MDP model
% Markov decision process is  1) observe state - 2) action (reserved) 3) Influence R 4) Dev 5) Inf Dev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
% G: adjacency matrix
% parentParcels: cellarray of parent ids (parentParcels{i,:} is the ids that
%                machine i depends on in the network)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%
% ddMDP: mdp struct (factored MDP) encoded with algebraic decision diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments:
%
% Generates a MDP struct (factored MDP) of an invasive species network problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% State variables %%%%%%%%%%%%%%%%%%%%%%%

%eff=tr.pm;
%noeff=1-tr.pm;
% detect0=tr.pdm;
% detect=tr.pds;
% nodetect=1-tr.pds;
%p=tr.pn;

%cost_m=1;
%reward=cost_m*RewM_ratio;
%cost_s=cost_m/MS_ratio;

ddMDP.nStateVars = length(parentParcels);

for i = 1:ddMDP.nStateVars
    ddMDP.stateVars(i).arity = 3;
    ddMDP.stateVars(i).name = sprintf('Site%i',i);
    ddMDP.stateVars(i).id = i;
    ddMDP.stateVars(i).valNames = {'avail','res', 'dev'};
end

%%%%%%%%%%%%%%%%%% Globals %%%%%%%%%%%%%%%%%%%%%%%

ddMDP.nVars = ddMDP.nStateVars;

Global.setVarDomSize([ddMDP.stateVars(:).arity, ddMDP.stateVars(:).arity]);

varNames = {ddMDP.stateVars(:).name};
for i = 1:ddMDP.nVars
    varNames = {varNames{:},[varNames{i},char(39)]};
end
Global.setVarNames(varNames);

for i = 1:ddMDP.nStateVars
    Global.setValNames(i,ddMDP.stateVars(i).valNames);
end

for i = 1:ddMDP.nStateVars
    Global.setValNames(ddMDP.nVars+i,ddMDP.stateVars(i).valNames);
end

% check here
%%%%%%%%%%%%%%%%% Transition function %%%%%%%%%%%%%%%%%%%%%
% manage parcel i
% do nothing
%
ddMDP.nActions = ddMDP.nStateVars + 1; % pas de SURVEY

for actId = 1:ddMDP.nActions
    if actId <= ddMDP.nStateVars
        ddMDP.actions(actId).name = sprintf('reserveS%i',actId);
    else
        ddMDP.actions(actId).name = 'do_Nothing';
    end

    for parcelId = 1:ddMDP.nStateVars
        % parcel is managed
        if actId == parcelId

            nr=0;   % number of reserved neighbours
            nd=0;   % number of developed neigbours
            liste = parentParcels{parcelId}; % copy the list of parents of node parcelId
            Left=BuildDDSIR1(liste,parcelId,pr,pd,pir,pid,nr,nd,0);  % p(site'=avail| ...
            Centre=BuildDDSIR1(liste,parcelId,pr,pd,pir,pid,nr,nd,1); % p(site'=reserv|...
            Right=BuildDDSIR1(liste,parcelId,pr,pd,pir,pid,nr,nd,2); % p(site'=dev|...


            ddMDP.actions(actId).transFn(parcelId)= DDnode.myNew(parcelId+ddMDP.nVars, ...
                [DDnode.myNew(parcelId,[Left,DDleaf.myNew(0),DDleaf.myNew(0)]),... % p(site'=avail| ...
                DDnode.myNew(parcelId,[Centre,DDleaf.myNew(1),DDleaf.myNew(0)]),...         % p(site'=reserv|...
                DDnode.myNew(parcelId,[Right,DDleaf.myNew(0),DDleaf.myNew(1)])]);           % p(site'=dev|...
        else % parcel not managed
            nr=0;   % number of reservd neighbours
            nd=0;   % number of developed neigbours
            prZ=0;
            liste = parentParcels{parcelId}; % copy the list of parents of node parcelId
            Left=BuildDDSIR1(liste,parcelId,prZ,pd,pir,pid,nr,nd,0);  % p(site'=avail| ...
            %Left.display()
            Centre=BuildDDSIR1(liste,parcelId,prZ,pd,pir,pid,nr,nd,1); % p(site'=reserv|...
            Right=BuildDDSIR1(liste,parcelId,prZ,pd,pir,pid,nr,nd,2); % p(site'=dev|...


            ddMDP.actions(actId).transFn(parcelId)= DDnode.myNew(parcelId+ddMDP.nVars, ...
                [DDnode.myNew(parcelId,[Left,DDleaf.myNew(0),DDleaf.myNew(0)]),... % p(site'=avail| ...
                DDnode.myNew(parcelId,[Centre,DDleaf.myNew(1),DDleaf.myNew(0)]),...         % p(site'=reserv|...
                DDnode.myNew(parcelId,[Right,DDleaf.myNew(0),DDleaf.myNew(1)])]);           % p(site'=dev|...
        end
    end
end

%%%%%%%%%%%%%%%%%%% Cost function %%%%%%%%%%%%%%%%%%%%%%
for actId = 1:ddMDP.nActions
    % cost for managing
    if actId <= ddMDP.nStateVars
        ddMDP.actions(actId).costFn=[DDleaf.myNew(cost)];
    else % action do nothing
        ddMDP.actions(actId).costFn=DD.zero;
    end
end

%%%%%%%%%%%%%%%%%%% Reward function %%%%%%%%%%%%%%%%%%%%%%
sitelist=[2:ddMDP.nStateVars];
Left=BuildDDRew(sitelist,BR,zeros(1,size(BR,2))); %available
Centre=BuildDDRew(sitelist,BR,BR(1,:)); %reserved
Right=BuildDDRew(sitelist,BR,zeros(1,size(BR,2))); %developed

ddMDP.rewFn=DDnode.myNew(1,[Left,Centre,Right]);

%%%%%%%%%%%%%%%%%% Discount Factor %%%%%%%%%%%%%%%%%%%%%%%%

ddMDP.discFact = 0.96;

%%%%%%%%%%%%%%%%%% Initial Belief State %%%%%%%%%%%%%%%%%%%%%%%%

ddInit = DD.one;
for parcelId = 1:ddMDP.nStateVars
    ddInit = OP.mult(ddInit, DDnode.myNew(parcelId,[DDleaf.myNew(1),DDleaf.myNew(0),DDleaf.myNew(0)]));
end
ddMDP.initialBelState = ddInit;

%%%%%%%%%%%%%%%%%% Tolerance %%%%%%%%%%%%%%%%%%%%%%%%

maxVal = -inf;
minVal = inf;
for actId = 1:ddMDP.nActions
    maxVal = max(maxVal,OP.maxAll(OP.addN(ddMDP.actions(actId).rewFn)));
    minVal = min(minVal,OP.minAll(OP.addN(ddMDP.actions(actId).rewFn)));
end
maxDiffRew = maxVal - minVal;
maxDiffVal = maxDiffRew/(1-min(0.95,ddMDP.discFact));
ddMDP.tolerance = 1e-5 * maxDiffVal;
%ddMDP.tolerance = 0.01;
%for actId = 1:ddMDP.nActions
%  ddMDP.actions(actId).rewFn = OP.approximate(ddMDP.actions(actId).rewFn, ddMDP.tolerance, [0]);
%end
