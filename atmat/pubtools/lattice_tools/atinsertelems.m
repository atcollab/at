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
% if ELEM = [], nothing is inserted, only the splitting takes place


if islogical(indef),refpts=find(refpts); end

deb=1;
slices=arrayfun(@insert,refpts,'UniformOutput',false);  %Split into sections
newring=cat(1,slices{:},ring(deb:end)); %Concatenate all sections

    function slice=insert(idx)
        slice=[ring(deb:idx-1);atsplitelem(ring{idx},varargin{:})];
        deb=idx+1;
    end
end
