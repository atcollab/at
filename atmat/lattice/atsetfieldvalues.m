function newring = atsetfieldvalues(ring,varargin)
%ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
%
%   Note that the calling syntax must be in the form of assignment:
%   RING = ATSETFIELDVALUES(RING,...)
%   MATLAB does not modify variables that only appear on the right
%   hand side as arguments.
%
%NEWRING=ATSETFIELDVALUES(RING,'field',VALUES)
%   In this mode, the function will set values on all the elements of RING
%
%NEWRING=ATSETFIELDVALUES(RING,INDEX,'field',VALUES)
%   In this mode, the function will set values on the elements of RING
%   specified by INDEX, given as a list of indices or as a logical mask
%
%NEWRING=ATSETFIELDVALUES(RING,'field',VALUESTRUCT)
%   In this mode, the function will set values on the elements of RING
%   whose family names are given by the field names of VALUESTRUCT
%
%NEWRING=ATSETFIELDVALUES(RING,RINGINDEX,....,VALUESTRUCT)
%   As in the previous mode, the function will set values on the elements
%   of RING whose family names are given by the field names of VALUESTRUCT.
%   But RINGINDEX=atindex(RING) is provided to avoid multiple computations.
%
% Field selection
% ---------------------------------------------------------
% NEWRING = ATSETFIELD(RING,'field',VALUES)
%   For each I=1:length(RING), set RING{I}.FIELD=value
%
% NEWRING = ATSETFIELD(RING,'field',{M,N},VALUES)
%   For each I=1:length(RING), set RING{I}.FIELD(M,N)=value
%
% More generally,
% NEWRING = ATSETFIELD(RING,subs1,subs2,...,VALUES)
%   For each I=1:length(RING), SETFIELD(RING{I},subs1,subs2,...,value)
%
% The first dimension of VALUES must be either length(INDEX) or 1 (the value
% will be repeated for each element). For a vector to be repeated, enclose
% it in a cell array.
%
% Value format
% ---------------------------------------------------------
% Cell array VALUES
% -----------------
% Mx1 or 1xM cell array : one cell per element
% 1x1 cell array : cell 1 is affected to all selected elements
%
% Character array VALUES
% ---------------------
% 1xN char array (string) : the string as affected to all selected elements
% MxN char array : one row per element
%
% Numeric array VALUES
% --------------------
% 1x1 (scalar) : the value is affected to all selected elements
% Mx1 (column) : one value per element
% MxN (matrix) : one row affected per element. If M==1, the single row
%                is affected to all selected elements
% To assign column vectors to parameters, use a cell array VALUES, with one
% column per cell
%
% See also ATGETFIELDVALUES ATGETCELLS SETCELLSTRUCT FINDCELLS

if isstruct(varargin{end})
    if isstruct(varargin{1})
        ringidx=varargin{1};
        args=varargin(2:end-1);
    else
        ringidx=atindex(ring);
        args=varargin(1:end-1);
    end
    newring=ring;
    for fn=fieldnames(varargin{end})'
        fname=char(fn);
        newring(ringidx.(fname))=atsetfield(newring(ringidx.(fname)),...
            args{:},varargin{end}.(fname));
    end
elseif islogical(varargin{1}) || isnumeric(varargin{1})
    newring=ring;
    newring(varargin{1})=atsetfield(newring(varargin{1}),varargin{2:end});
else
    newring=atsetfield(ring,varargin{:});
end

    function newring = atsetfield(ring,varargin)        
        if ischar(varargin{end})
            values=cellstr(varargin{end});
        elseif isnumeric(varargin{end})
            values=num2cell(varargin{end},2);
        else
            values=varargin{end};
        end
        if isscalar(values)
            values=values(ones(length(ring),1));
        end
        newring=cellfun(@(el,value) setfield(el,varargin{1:end-1},squeeze(value)),...
            ring,values(:),'UniformOutput',false); %#ok<SFLD>
    end

end
