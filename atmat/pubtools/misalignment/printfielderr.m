function varargout = printfielderr(varargin)
% PRINTFIELDERR will give a short display of the field errors being
% applied to the ring.
%
% PRINTFIELDERR(['famname', 'ind','file',FILENAME]) where 'ind' will print out the
% field error for each element and the 'file' option and FILENAME have
% to be specified together in order to output the data to a file. 'famname'
% is just a list of family elements that you're interested in and has to
% come before 'ind' or 'file'.
%
% Eg. printfielderr q1 q2 q3 ind file ferr.txt;
%     printfielderr('q1','q2','q3','ind','file','ferr.txt');
%     printfielderr({'q1','q2','q3'},'ind','file','ferr.txt');

global THERING FAMLIST

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

ferr = getappdata(0,'FielderrData');

if isempty(ferr)
    disp('No field error data found. See SETFIELDERR for more info');
    return
elseif ferr.currseed == 0
    disp('Field error data that calculated but not applied yet. See APPLYFIELDERR.');
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
    fprintf(fid,'Printing individual magnet field errors\n');
    fprintf(fid,'Current seed being used = %d\n',ferr.currseed);
    fprintf(fid,'========================================\n');
    fprintf(fid,' Element name            |  Multipole errors     |\n');
    for i=1:length(ferr.data)
        if ~isempty(strmatch(ferr.data(i).name,strvcat(ellist)))
            if isempty(ferr.data(i).valA)
                fprintf(fid,'%14s\n',ferr.data(i).name);
            else
                indA = find(ferr.data(i).valA(:,ferr.currseed) ~= 0);
                indB = find(ferr.data(i).valB(:,ferr.currseed) ~= 0);
                fprintf(fid,'%14s (PolynomA): ',ferr.data(i).name);
                if ~isempty(indA)
                    fprintf(fid,'(pole%02d) %+.3e  |',[indA; ferr.data(i).valA(indA,ferr.currseed)]);
                end
                fprintf(fid,'\n%14s (PolynomB): ',' ');
                if ~isempty(indB)
                    fprintf(fid,'(pole%02d) %+5.3e  |',[indB; ferr.data(i).valB(indB,ferr.currseed)]);
                end
                fprintf(fid,'\n');
            end
        end
    end
else
    % summary groups by family
    fprintf(fid,'Printing SUMMARY of the magnet misalignments\n');
    fprintf(fid,'Current seed being used = %d\n',ferr.currseed);
    fprintf(fid,'========================================\n');
    fprintf(fid,'Element name  | pole  | mean (PolynomA, PolynomB)|  std (PolynomA, PolynomB)\n');
    for i=1:nel
        summary(i).name = ellist{i};
        indices = findcells(THERING,'FamName',ellist{i});
        if ~isempty(indices) & ferr.data(indices(1)).used
            % Cancatenate all misalignments into one large matrix
            cat_ferr_vectorsA = cell2mat({ferr.data(indices).valA});
            cat_ferr_vectorsB = cell2mat({ferr.data(indices).valB});
            % Select only the field error vector for current seed
            ferr_vec_for_currentseedA = cat_ferr_vectorsA(:,ferr.currseed:ferr.numseeds:end);
            ferr_vec_for_currentseedB = cat_ferr_vectorsB(:,ferr.currseed:ferr.numseeds:end);
            % Calculate mean along rows
            summary(i).meanA = mean(ferr_vec_for_currentseedA,2);
            summary(i).deviationA = std(ferr_vec_for_currentseedA,0,2);
            summary(i).meanB = mean(ferr_vec_for_currentseedB,2);
            summary(i).deviationB = std(ferr_vec_for_currentseedB,0,2);
            
            fprintf(fid,'%12s', summary(i).name);
            fprintf('   pole%02d : %+5.3e, %+5.3e : %+5.3e, %+5.3e\n',1,...
                summary(i).meanA(1),summary(i).meanB(1),...
                summary(i).deviationA(1),summary(i).deviationB(1));
            
            indA = [2:length(summary(i).meanA)];
            indB = [2:length(summary(i).meanB)];
            fprintf('               pole%02d : %+5.3e, %+5.3e : %+5.3e, %+5.3e\n',...
                [indA; summary(i).meanA(indA)'; summary(i).meanB(indB)';...
                summary(i).deviationA(indA)';summary(i).deviationB(indB)']);
        else
            summary(i).mean = [];
            summary(i).deviation = [];
        end
    end

end

if fid > 2
    fclose(fid);
end