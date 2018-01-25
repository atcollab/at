function values = getcellstruct(CELLARRAY,field,index,varargin)
%GETCELLSTRUCT retrieves the field values MATLAB cell array of structures
%
% VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX,M,N)
%
% VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX,M) can be used
%   for one dimensional vectors
%
% VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX) is the same as
%   GETCELLSTRUCT(CELLARRAY,'field',index,1,1) if the field data
%   is a scalar
%
% VALUES = GETCELLSTRUCT(CELLARRAY,'field',INDEX) is a MATLAB cell array
% 	 of strings if specified fields contain strings.
%
% See also ATGETFIELDVALUES SETCELLSTRUCT FINDCELLS

if isempty(varargin)
    values=atgetfieldvalues(CELLARRAY(index),field);
elseif length(varargin)==1
    values=atgetfieldvalues(CELLARRAY(index),field,varargin{1});
else
    values=atgetfieldvalues(CELLARRAY(index),field,{varargin{1},varargin{2}});
end
