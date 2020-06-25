function [graphM,terminalAct]=consP_build_policy_tree(n,filename)
% input:
% n is the number of node
% filename is a policy dot file from spudd
% output:
% graphM is a (m,4) matrix representing the policy tree
% for terminal node, actions can be found in terminalAct

[y,Cnode,terminalAct, sitesID]=find_node_label(filename,'{ rank');
[yg,Clink]=find_graph_label(filename,' -> ');

[poub n]=size(Cnode);
[poub m]=size(Clink);

graphM=zeros(n,4);

% go trough each link and populate the matrix
for i=1:m
    inode=str2num(char(Clink{i}(1)))+1;
    next_node=str2num(char(Clink{i}(2)))+1;
    str_link=char(Clink{i}(3));

    if strcmp(str_link,'avail')
        graphM(inode,1)= next_node;
    elseif strcmp(str_link,'res')
        graphM(inode,2)= next_node;
    elseif strcmp(str_link,'dev')
        graphM(inode,3)= next_node;
    end
end
% we have the links we need the actions
graphM(:,4)=sitesID(:,2);
end

function [y,Ct, terminal, sites] = find_node_label(filename, literal)
% Search for number of string matches per line.
% return y the number of nodes of the graph
% and Ct their label for example Ct{4} contains '3'    'reserveS2'
% terminal yes|no indicate if a node is terminal

fid = fopen(filename);
if fid==-1
    error('file does not exist')
else
y = 0;
tline = fgetl(fid);
while ischar(tline)
   matches = strfind(tline, literal);

   num = length(matches);
   if num > 0
       expr=';".(\d*)" [label\s?=\s?"(\w*)\W?.*?\s?"';
       [tok mat] = regexp(tline, expr, 'tokens', 'match');
       y = y + num;
       %fprintf(1,'%d:%s\n',num,tline);
       if ~isempty(tok)
        Ct{y}=tok{1,:};

        % action nodes
        x=char(Ct{y}(2));
        [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(x,'reserveS');
        if ~isempty(start_idx)
            expr='reserveS(\d*)';       % in case nb sites >9
            [tok mat] = regexp(x, expr, 'tokens', 'match');
            terminal(y)=str2num(char(tok{:}));
            %terminal(y)=str2num(x(end_idx+1));
        else
            [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(x,'do_Nothing');
            if ~isempty(start_idx)
                terminal(y)=0;  % do nothing
            else
                terminal(y)=-1; % not terminal -> keep going
            end
        end

        % site nodes ID
       expr='Site(\d*)';
       [tok mat] = regexp(tline, expr, 'tokens', 'match');
       [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(x,'Site');
        if ~isempty(start_idx)
            sites(y,1)=str2num(char(Ct{y}(1)));
            sites(y,2)=str2num(char(tok{:}));
        else
            sites(y,1)=str2num(char(Ct{y}(1)));
            sites(y,2)=-1;
        end

       else
           'error<<'
       end
   end
   tline = fgetl(fid);
end
fclose(fid);
end
end

function [y,Ct] = find_graph_label(filename, literal)
% Search for number of string matches per line.
% return y the number of links and their attribute
% for example Ct{1} contains '2'    '3'    'avail'

fid = fopen(filename);
y = 0;
tline = fgetl(fid);
while ischar(tline)
   matches = strfind(tline, literal);

   num = length(matches);
   if num > 0
       expr='".(\d*)"\s->\s".(\d*)"\s[label\s=\s"(\w*)"';
       [tok mat] = regexp(tline, expr, 'tokens', 'match');
       y = y + num;

      %fprintf(1,'%d:%s\n',num,tline);
       if ~isempty(tok)
        Ct{y}=tok{1,:};
       else
           'error<<'
       end
   end
   tline = fgetl(fid);
end
fclose(fid);
end
