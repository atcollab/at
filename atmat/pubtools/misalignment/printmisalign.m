function varargout = printmisalign(varargin)
% PRINTMISALIGN will give a short display of the misalignements being
% applied to the ring.
%
% PRINTMISALIGN(['famname', 'ind','file',FILENAME]) where 'ind' will pring out the
% misalignment for each element and the 'file' option and FILENAME have
% to be specified together in order to output the data to a file. 'famname'
% is just a list of family elements that you're interested in and has to
% come before 'ind' or 'file'.
%
% Eg. printmisalign q1 q2 q3 ind file mis.txt;
%     printmisalign('q1','q2','q3','ind','file','mis.txt');
%     printmisalign({'q1','q2','q3'},'ind','file','mis.txt');

global THERING FAMLIST

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
elseif mis.currseed == 0
    disp('Alignment data thats calculated not applied yet. See APPLYCMISALIGN.');
    return
end

ind = 0; output = 0; filename = '';
ellist = {}; nel = 0;
for i=1:nargin
    if ischar(varargin{i})
        switch varargin{i}
            case 'ind'
                ind = 1;
            case 'file'
                output = 1;
                filename = varargin{i+1};
                varargin{i+1} = 0;
            otherwise
                nel = nel + 1;
                ellist{nel} = varargin{i};
        end
    elseif i == 1 & iscell(varargin{i})
        nel = length(varargin{i});
        ellist = varargin{i};
    end
end

% If no elements are specified then all elements will be printed out
if nel == 0
    for i=1:length(FAMLIST)
        nel = nel + 1;
        ellist{nel} = FAMLIST{i}.FamName;
    end
end

if output
    fid = fopen(filename,'w');
    if fid == 0
        error(['Could not open file: ' filename]);
        return
    end
else
    % Print to screen
    fid = 1;
end

fprintf(fid,'%s\n',datestr(clock,31));
if ind
    % individual elements
    fprintf(fid,'Printing individual magnet misalignments\n');
    fprintf(fid,'Current seed being used = %d\n',mis.currseed);
    fprintf(fid,'========================================\n');
    fprintf(fid,' Element name    |  x  |  rx  |  y  |  ry  |  s  |  rs \n');
    for i=1:length(mis.data)
        if ~isempty(strmatch(mis.data(i).name,strvcat(ellist)))
            if isempty(mis.data(i).val)
                fprintf(fid,'%14s\n',mis.data(i).name);
            else
                fprintf(fid,'%14s | %6.3e  | %6.3e | %6.3e | %6.3e | %6.3e | %6.3e\n',...
                    mis.data(i).name, mis.data(i).val(1,mis.currseed), ...
                    mis.data(i).val(2,mis.currseed), ...
                    mis.data(i).val(3,mis.currseed), ...
                    mis.data(i).val(4,mis.currseed), ...
                    mis.data(i).val(5,mis.currseed), ...
                    mis.data(i).val(6,mis.currseed));
            end
        end
    end
else
    % summary groups by family
    fprintf(fid,'Printing SUMMARY of the magnet misalignments\n');
    fprintf(fid,'Current seed being used = %d\n',mis.currseed);
    fprintf(fid,'========================================\n');
    fprintf(fid,'Element name  | coord |    mean    |  deviation  \n');
    for i=1:nel
        summary(i).name = ellist{i};
        indices = findcells(THERING,'FamName',ellist{i});
        if ~isempty(indices) & mis.data(indices(1)).used
            % Cancatenate all misalignments into one large matrix
            cat_mis_vectors = cell2mat({mis.data(indices).val});
            % Select only the misalignment vector for current seed
            mis_vec_for_currentseed = cat_mis_vectors(:,mis.currseed:mis.numseeds:end);
            % Calculate mea along rows
            summary(i).mean = mean(mis_vec_for_currentseed,2);
            summary(i).deviation = std(mis_vec_for_currentseed,0,2);

            fprintf(fid,'%10s        x   : %6.3e : %6.3e\n',summary(i).name,...
                summary(i).mean(1),summary(i).deviation(1));
            fprintf(fid,'                 rx   : %6.3e : %6.3e\n',...
                summary(i).mean(2),summary(i).deviation(2));
            fprintf(fid,'                  y   : %6.3e : %6.3e\n',...
                summary(i).mean(3),summary(i).deviation(3));
            fprintf(fid,'                 ry   : %6.3e : %6.3e\n',...
                summary(i).mean(4),summary(i).deviation(4));
            fprintf(fid,'                  s   : %6.3e : %6.3e\n',...
                summary(i).mean(5),summary(i).deviation(5));
            fprintf(fid,'                 rs   : %6.3e : %6.3e\n',...
                summary(i).mean(6),summary(i).deviation(6));
        else
            summary(i).mean = [];
            summary(i).deviation = [];
        end
    end

end

if fid > 2
    fclose(fid);
end