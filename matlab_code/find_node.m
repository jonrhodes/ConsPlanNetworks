function [terminal,sites] = find_node(filename)
% Search for number of string matches per line.
% terminal(x,y)=[nodeid, {-1 if node is not terminal,0=do_nothing, action}] indicate if a node is terminal, and what
% action to perform
% sites(x,y)=[nodeid, {siteid, -1 if node is an action}]

fid = fopen(filename);
if fid==-1
    error('file does not exist')
else
    y = 0;
    tline = fgetl(fid);
    while ischar(tline)
        % site nodes
        matches = strfind(tline, '[shape=ellipse');  % site nodes
        num = length(matches);
        if num > 0   % we have a match!
            expr=';".(\d*)" [label="Site(\d*)';
            [tok mat] = regexp(tline, expr, 'tokens', 'match');
            y = y + 1;
            %fprintf(1,'%d:%s\n',y,tline);
            Ct{y}=tok{1,:};
            sites(y,1)=str2num(char(Ct{y}(1)));
            sites(y,2)=str2num(char(Ct{y}(2)));
            terminal(y,1)=sites(y,1); % not terminal -> keep going
            terminal(y,2)=-1; % not terminal -> keep going
        end

        %action nodes can be reserve site or no nothing
        matches = strfind(tline, '[shape=box');  %action nodes
        num = length(matches);
        if num > 0
            %fprintf(1,'%d:%s\n',y,tline);
            y = y + 1;
            expr=';".(\d*)" [label = "reserveS(\d*)';
            [tok mat] = regexp(tline, expr, 'tokens', 'match');

            if ~isempty(tok)
                Ct{y}=tok{1,:};
                terminal(y,1)=str2num(char(Ct{y}(1)));
                terminal(y,2)=str2num(char(Ct{y}(2)));
                sites(y,1)=terminal(y,1);
                sites(y,2)=-1;
            end
            expr=';".(\d*)" [label = "do_Nothing';
            [tok mat] = regexp(tline, expr, 'tokens', 'match');

            if ~isempty(tok)
                Ct{y}=tok{1,:};
                terminal(y,1)=str2num(char(Ct{y}(1)));  % do nothing
                terminal(y,2)=0;  % do nothing
                sites(y,1)=terminal(y,1);
                sites(y,2)=-1;
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

end
