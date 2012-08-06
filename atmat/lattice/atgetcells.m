function ok=atgetcells(cellarray, field, varargin)
%FINDCELLS performs a search on MATLAB cell arrays of structures
%
% INDEX = ATGETCELLS(RING, 'field')
%   returns indexes of elements that have a field named 'field'
%
% INDEX = ATGETCELLS(RING, 'field', VALUE1...)
%   returns indexes of elements whose field 'field'
%   is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
%   character strings or a number. If its a character string REGULAR
%   expressions can be used.
%
% INDEX = ATGETCELLS(RING, 'field', @TESTFUNCTION,...)
%   Uses the user-defined TESTFUNCTION to select array elements
%   TESTFUNCTION must be of the form:
%       OK=TESTFUNTION(ATELEM,FIELD,ARGS...)
%           OK : field FIELD of element ATELEM matches the condition
%               defined by ARGS
%
% OK is a logical array with the same size as RING, refering to matching
% elements in RING
%
% See also GETCELLSTRUCT, SETCELLSTRUCT, REGEXPI

% Check if the first argument is the cell array of structures
if(~iscell(cellarray) || ~isstruct(cellarray{1}) || isempty(cellarray))
    error('The first argument must be a non-empty cell array of structures')
end
% Check if the second argument is a string
if(~ischar(field))
    error('The second argument must be a character string')
end

NE = numel(cellarray);
ok=false(size(cellarray));

if nargin<3
    tesfunc=@(elem,field) true;
    vals={};
elseif isa(varargin{1},'function_handle')
    tesfunc=varargin{1};
    vals=varargin(2:end);
else
    tesfunc=@defaultfunc;
    vals=varargin;
end

for elem = 1:NE
    if isfield(cellarray{elem},field)
        ok(elem)=tesfunc(cellarray{elem},field,vals{:});
    end
end

    function ok=defaultfunc(el,field,varargin)
        ok=false;
        for j=1:length(varargin)
            if ischar(el.(field)) && ischar(varargin{j})
                if ~isempty(regexpi(el.(field),['^' varargin{j} '$']))
                    ok=true;
                    break;
                end
            elseif  isnumeric(el.(field)) && isnumeric(varargin{j})
                if el.(field) == varargin{j}
                    ok=true;
                    break;
                end
            end
        end
    end
end
