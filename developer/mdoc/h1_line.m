function hh = h1_line(filename)
%H1_LINE get the H1 line for a file

[~,name,~] = fileparts(filename);
h1 = ''; % default output
fid = fopen(filename); % open file
if fid ~= -1
    tline = fgetl(fid); % read first line
    while ischar(tline)
        k = strfind(tline,'%'); % find comment
        if ~isempty(k) % if it is found
            k = k(1);
            ispercents = false(size(tline(k:end)));
            ispercents(strfind(tline(k:end),'%'))=true;
            start = k+find(~(isspace(tline(k:end)) | ispercents),1,'first')-1;
            if ~isempty(start)
                tline = tline(start:end); % remove leading space/percent
                IX = strfind(lower(tline),lower(name));
                if isempty(IX)
                    error('H1:Invalid','%s: Invalid H1 line',name);
                else
                    if IX(1)==1
                        tline = tline(length(name)+1:end); % remove function name
                    end
                    tline = strtrim(tline); % remove any leading/trailing space
                end
                h1 = tline;
                h1 = strtrim(h1);
                if ~isempty(h1)
                    if strcmp(h1(end),'.') % remove trailing period
                        h1 = h1(1:end-1);
                    end
                    h1(1) = upper(h1(1)); % capitalize first letter
                end
            end
            tline = -1; % set tline to numeric
        else
            tline = fgetl(fid); % read next line
        end
    end
    fclose(fid);
    hh=struct('name',name,'h1',h1);
else
    error('AT:NoSuchFile','Cannot open file %s',filename);
end
end
