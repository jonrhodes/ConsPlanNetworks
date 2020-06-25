function action=getAction(Mstate,Mgraph,terminalAct)
% input:
% Mstate is a vector of size the number of sites and takes {1,2,3} value
% for each element. assuming: avail=1; res=2; dev = 3;
% e.g. state=[1 2 3] is site1 = avail, site2=res, site3 =dev
% output:
% action takes value [0,max number of sites] 0 is do nothing and x>0 is the
% site to reserve
% for example if action=10 => reserve site 10

% we have to start with the site at the root of the Mgraph

[site_r,node_r]=get_root(Mgraph);
k=length(Mstate);

for i=1:k % go through each state unless node is terminal
    if terminalAct(node_r,2)==-1  % node is a site
        val=Mstate(site_r); % avail=1, res=2, dev=3
        next_node=Mgraph(node_r,val);
        next_site=Mgraph(next_node,4);
        node_r=next_node;
        site_r=next_site;
    else %node is terminal => get action!
        break
    end

end
action=terminalAct(node_r,2); % in case action is not allocated prior
                            % to knowing all the states - e.g. normal case

if action>0
   if Mstate(action)~=1
       'bad action'
   end
end

if action<0
    'err getAction'
    disp(node_r)
    disp(Mstate)
end

end
