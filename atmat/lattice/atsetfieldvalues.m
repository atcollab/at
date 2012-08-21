function ring = atsetfieldvalues(ring,index,varargin)
%ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
%
%   Note that the calling syntax must be in the form of assignment:
%   RING = ATSETFIELDVALUES(RING,...)
%   MATLAB does not modify variables that only appear on the right
%   hand side as arguments.
%
% Numeric data
% ---------------------------------------------------------
% NEWRING = ATSETFIELDVALUES(RING,INDEX,'field',VALUES)
%   For each I=INDEX, set RING{INDEX(I)}.FIELD=VALUES(:,:,I)
%
% NEWRING = ATSETFIELDVALUES(RING,INDEX,'field',{M,N},VALUES)
%   For each I=INDEX, set RING{INDEX(I)}.FIELD(M,N)=VALUES(I)
%
% More generally,
% NEWRING = ATSETFIELDVALUES(RING,INDEX,subs1,subs2,...,VALUES)
%   For each I=INDEX, SETFIELD(RING{INDEX(I)},subs1,subs2,...,VALUES(...,I))
%
% The last dimension of VALUES must be either length(INDEX) or 1 (the value
% will be repeated for each element). For a vector to be repeated, enclose
% it in a cell array.
%
% Character array
% --------------------------------------------------------------------
% NEWRING = ATSETFIELDVALUES(RING,INDEX,'field',VALUES)
%   For each I=INDEX, set RING{INDEX(I)}.FIELD=VALUES(I,:)
%
% Each row of VALUES is affected to the specified field of each selected
% element of RING. The number of rows on values must be either
% length(INDEX) or 1 (the value will be repeated for each element).
%
% Cell array
% --------------------------------------------------------------------
% NEWRING = ATSETFIELDVALUES(RING,INDEX,'field',VALUES)
%   For each I=INDEX, set RING{INDEX(I)}.FIELD=VALUES{I}
%
% Each cell of VALUES is affected to the specified field of each selected
% element of RING.The length of VALUES must be either length(INDEX) or 1
% (the value will be repeated for each element). 
%
% See also ATGETFIELDVALUES ATGETCELLS SETCELLSTRUCT FINDCELLS

if ~iscell(ring) || isempty(ring) || ~isstruct(ring{1})
    error('The first argument must be a non-empty cell array of structures')
end

if islogical(index)
    index=find(index);
end

if ischar(varargin{end})            % Character array
    values=cellstr(varargin{end});
elseif isnumeric(varargin{end})     % Numeric data
    [r c]=size(varargin{end}); %#ok<ASGLU>
    if c==1
        values=num2cell(varargin{end});
    else
        values=squeeze(num2cell(varargin{end},1:ndims(varargin{end})-1));
    end
else                                % Cell array of strings
    values=varargin{end};
end

nvals = numel(index);
if length(values)==1
    values=values(ones(nvals,1));
end
for idv=1:nvals
    idx=index(idv);
    ring{idx}=setfield(ring{idx},varargin{1:end-1},values{idv}); %#ok<SFLD>
end
