function values = atgetfieldvalues(ring,varargin)
%ATGETFIELDVALUES retrieves the field values AT cell array of elements
%
% VALUES = ATGETFIELDVALUES(RING,'field') extracts the values of
% the field 'field' in all the elements of RING
%
% VALUES = ATGETFIELDVALUES(RING,INDEX,'field') extracts the values of
% the field 'field' in the elements of RING selected by INDEX
%
% if RING{I}.FIELD is a numeric scalar
%    VALUES is a length(INDEX) x 1 array
% otherwise
%    VALUES is a length(INDEX) x 1 cell array
%
%
% More generally ATGETFIELDVALUES(RING,INDEX,subs1,subs2,...) will call
%  GETFIELD(RING{I},subs1,subs2,...) for I in INDEX
%
% Examples:
%
% V=ATGETFIELDVALUES(RING,1:10,'PolynomB') is a 10x1  cell array
% such that V{I}=RING{I}.PolynomB for I=1:10
%
% V=ATGETFIELDVALUES(RING(1:10),'PolynomB',{1,2}) is a 10x1 array
% such that V(I)=RING{I},PolynomB(1,2)
%
%
% See also ATSETFIELDVALUES ATGETCELLS GETCELLSTRUCT FINDCELLS

if islogical(varargin{1}) || isnumeric(varargin{1})
    values=atgetfield(ring(varargin{1}),varargin{2:end});
else
    values=atgetfield(ring,varargin{:});
end

    function values = atgetfield(line,varargin)
        [values,isnumscal,isok]=cellfun(@scan,line(:),'UniformOutput',false);
        isok=cat(1,isok{:});
        if all(cat(1,isnumscal{isok}))
            values(~isok)={NaN};
            values=cat(1,values{:});
        end
        
        function [val,isnumscal,isok]=scan(el)
            try
                val=getfield(el,varargin{:});
                isnumscal=isnumeric(val) && isscalar(val);
                isok=true;
            catch
                val=[];
                isnumscal=false;
                isok=false;
            end
        end
    end

end
