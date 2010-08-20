function varargout = undofielderr(varargin)
% UNDOMISALIGN will 'zero' all field errors in THERING 

global THERING

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

[reg prop] = parseparams(varargin);
verbatim = 1;
for i=1:length(prop)
    if strcmpi(prop,'quiet')
        verbatim = 0;
    end
end

ferr = getappdata(0,'FielderrData');

if isempty(ferr)
    disp('No field error data found. See SETFIELDERR for more info');
    return
end


if isfield(ferr,'data') & ferr.currseed > 0
    % Indices of elements where misalignemnt data may need to be applied.
    indices = find(cell2mat({ferr.data.used}));
    for i=indices
        THERING{i}.PolynomA = THERING{i}.PolynomA - ferr.data(i).valA(:,ferr.currseed)';
        THERING{i}.PolynomB = THERING{i}.PolynomB - ferr.data(i).valB(:,ferr.currseed)';
        if isfield(THERING{i},'K')
            THERING{i}.K = THERING{i}.PolynomB(2);
        end
    end
end

if verbatim
    disp('Field errors have been reset. Please still check.');
end

