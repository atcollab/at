function newring=atinsertelems(ring,refpts,varargin)
%ATINSERTELEMS Insert elements at given locations in a line
%
%NEWLINE=ATINSERTELEMS(LINE,REFPTS,FRAC1,ELEM1[,FRAC2,ELEM2...])
%   a new line is created by inserting elements at each location specified
%   by REFPTS.
%
%LINE:      Cell array of structures
%REFPTS:    Insertion points (index list or logical mask)
%FRAC:      Location of the inserted element ELEM within LINE{REFPTS(i)}
%           0<=FRAC<=1
% if FRAC = 0, ELEM is inserted before LINE{REFPTS(i)} (no splitting)
% if FRAC = 1, ELEM is inserted after LINE{REFPTS(i)} (no splitting)
% if FRAC = NaN, LINE{REFPTS(i)} is replaced by ELEM (no check for identical length)
% if ELEM = [], nothing is inserted, only the splitting takes place
%
% FRAC and ELEM must be scalars or array of the same size as REFPTS
%
% See also ATSPLITELEM ATDIVELEM

if islogical(refpts),refpts=find(refpts); end

vargs=cellfun(@unfold,varargin,'UniformOutput',false);

deb=1;
slices=cellfun(@insert,num2cell(refpts),vargs{:},'UniformOutput',false);  %Split into sections
newring=cat(1,slices{:},ring(deb:end)); %Concatenate all sections

    function slice=insert(idx,varargin)
        if isnan(varargin{1})
            if iscell(varargin{2})
                slice=[ring(deb:idx-1);varargin{2}];
            else
                slice=[ring(deb:idx-1);varargin(2)];
            end
        else
            slice=[ring(deb:idx-1);atsplitelem(ring{idx},varargin{:})];
        end
        deb=idx+1;
    end

    function val=unfold(arg)
        if isempty(arg)
            val=cell(size(refpts));
        elseif isscalar(arg)
            val=arg(ones(size(refpts)));
        else
            val=reshape(arg,size(refpts));
        end
        if ~iscell(arg)
            val=num2cell(val);
        end
    end
end
