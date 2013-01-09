function newring = atsetfieldvalues(ring,varargin)
%ATSETFIELDVALUES sets the field values of MATLAB cell array of structures
%
%   Note that the calling syntax must be in the form of assignment:
%   RING = ATSETFIELDVALUES(RING,...)
%   MATLAB does not modify variables that only appear on the right
%   hand side as arguments.
%
%NEWRING=ATSETFIELDVALUES(RING,INDEX,....,VALUES)
%   In this mode, the function will set values on the elements of RING
%   specified by INDEX, given as a list of indices or as a logical mask
%
%NEWRING=ATSETFIELDVALUES(RING,....,VALUESTRUCT)
%   In this mode, the function will set values on the elements of RING
%   whose family names are given by the field names of VALUESTRUCT
%
%NEWRING=ATSETFIELDVALUES(RING,RINGINDEX,....,VALUESTRUCT)
%   As in the previous mode, the function will set values on the elements
%   of RING whose family names are given by the field names of VALUESTRUCT.
%   But RINGINDEX is provided to avoid multiple computation.
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

if islogical(varargin{1})
    newring=doset(ring,find(varargin{1}),varargin{2:end});
elseif isnumeric(varargin{1})
    newring=doset(ring,varargin{:});
elseif isstruct(varargin{end})
    if isstruct(varargin{1})
        ringidx=varargin{1};
        ib=2;
    else
        ringidx=atindex(ring);
        ib=1;
    end
    newring=ring;
    for fn=fieldnames(varargin{end})'
        fname=char(fn);
        newring=doset(newring,ringidx.(fname),varargin{ib:end-1},varargin{end}.(fname));
    end
end

    function ring=doset(ring,index,varargin)
        if ischar(varargin{end})            % Character array
            values=cellstr(varargin{end});
        elseif isnumeric(varargin{end})     % Numeric data
            [r,c]=size(varargin{end}); %#ok<ASGLU>
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
    end

end
