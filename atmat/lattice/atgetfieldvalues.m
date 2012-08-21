function values = atgetfieldvalues(ring,index,varargin)
%ATGETFIELDVALUES retrieves the field values AT cell array of elements
%
% VALUES = ATGETFIELDVALUES(RING,INDEX,'field') extracts the values of
% the field 'field' in the elements of RING specified by INDEX
%
% IF RING{I}.FIELD is a string, VALUES is a length(INDEX) x 1 cell array of
% strings
%
% IF RING{I}.FIELD is numeric, VALUES is a multi-dimensional array such
% that values(:,:,I) is equal to RING{I}.FIELD
%
% More generally ATGETFIELDVALUES(RING,INDEX,subs1,subs2,...) will call
%  GETFIELD(RING{I},subs1,subs2,...) for I in INDEX
%
% Example:
% V=ATGETFIELDVALUES(RING,1:10,'PolynomB') is a 1xNx10 array
% such that V(:,:,I)=RING{INDEX(I)}.PolynomB for I=1:10
%
% V=ATGETFIELDVALUES(RING,1:10,'PolynomB',{1,2}) is a 1x1x10 array
% such that V(1,1,I)=RING{INDEX(I)},PolynomB(1,2)
%
%
% See also ATSETFIELDVALUES ATGETCELLS GETCELLSTRUCT FINDCELLS

if ~iscell(ring) || isempty(ring) || ~isstruct(ring{1})
    error('The first argument must be a non-empty cell array of structures')
end

if islogical(index)
    index=find(index);
end

nvals = numel(index);

vals=cell(nvals,1);
vnum=true;
%vchar=true;
dims=[0 0];
for idv=1:nvals
    try
        vv=getfield(ring{index(idv)},varargin{:});
        vals{idv}=vv;
        if isnumeric(vv)
%            vchar=false;
            dims=max([dims;size(vv)]);
            classnam=class(vv);
        else
            vnum=false;
        end
        
    catch %#ok<*CTCH>
    end
end

if vnum
    values = NaN([dims nvals],classnam);
    for idv = 1:nvals
        try
            [vl,vc]=size(vals{idv});
            values(1:vl,1:vc,idv) = vals{idv};
        catch
        end
    end
else
    values=vals;
end


