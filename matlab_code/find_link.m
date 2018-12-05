function [y,Ct] = find_link(filename)
% Search for number of string matches per line.
% return y the number of links and their attribute
% for example Ct{1} contains '2'    '3'    'avail'

fid = fopen(filename);
y = 0;
tline = fgetl(fid);
while ischar(tline)
    matches = strfind(tline, '->');
    num = length(matches);
    if num > 0
        expr='"a(\d*)"\s->\s"a(\d*)"\s[label\s=\s"(\w*)"';
        [tok mat] = regexp(tline, expr, 'tokens', 'match');
        y = y + 1;
        %fprintf(1,'%d:%s\n',y,tline);
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