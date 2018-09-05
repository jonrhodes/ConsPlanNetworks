function ddMDP=fillddMDPNetSIR3(pr,pd,pir,pid,cost,BR,parentParcels)

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
        % parcel is managed => Look at the neighbors too!!!
        if actId == parcelId
            nr=0;   % number of reserved neighbours
            nd=0;   % number of developed neigbours
            liste = parentParcels{parcelId}; % copy the list of parents of node parcelId
            Left=BuildDDSIR3(liste,parcelId,pr,pd,pir,pid,nr,nd,0);  % p(site'=avail| ...
            % check BUG pir and nr, nd inverted in other file.
            Centre=BuildDDSIR3(liste,parcelId,pr,pd,pir,pid,nr,nd,1); % p(site'=reserv|...
            Right=BuildDDSIR3(liste,parcelId,pr,pd,pir,pid,nr,nd,2); % p(site'=dev|...
            
            ddMDP.actions(actId).transFn(parcelId)= DDnode.myNew(parcelId+ddMDP.nVars, ...
                [DDnode.myNew(parcelId,[Left,DDleaf.myNew(0),DDleaf.myNew(0)]),... % p(site'=avail| ...
                DDnode.myNew(parcelId,[Centre,DDleaf.myNew(1),DDleaf.myNew(0)]),...         % p(site'=reserv|...
                DDnode.myNew(parcelId,[Right,DDleaf.myNew(0),DDleaf.myNew(1)])]);           % p(site'=dev|...
            
        else % parcel not managed
            nr=0;   % number of reserved neighbours
            nd=0;   % number of developed neigbours
            prZ=0;
            pirZ=0;
            liste = parentParcels{parcelId}; % copy the list of parents of node parcelId
            Left=BuildDDSIR3(liste,parcelId,prZ,pd,pirZ,pid,nr,nd,0);  % p(site'=avail| ...
            % check BUG pir and nr, nd inverted in other file.
            Centre=BuildDDSIR3(liste,parcelId,prZ,pd,pirZ,pid,nr,nd,1); % p(site'=reserv|...
            Right=BuildDDSIR3(liste,parcelId,prZ,pd,pirZ,pid,nr,nd,2); % p(site'=dev|...
            
            ddMDP.actions(actId).transFn(parcelId)=DDnode.myNew(parcelId+ddMDP.nVars, ...
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
    maxVal = max(maxVal,OP.maxAll(OP.addN(ddMDP.actions(actId).costFn)));
    minVal = min(minVal,OP.minAll(OP.addN(ddMDP.actions(actId).costFn)));
end
maxDiffRew = maxVal - minVal;
maxDiffVal = maxDiffRew/(1-min(0.95,ddMDP.discFact));
ddMDP.tolerance = 1e-5 * maxDiffVal;

