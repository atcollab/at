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


% See also SETCELLSTRUCT FINDCELLS 
if(~iscell(CELLARRAY) | ~isstruct(CELLARRAY{1}) | isempty(CELLARRAY))
   error('The first argument must be a non-empty cell array of structures') 
end
% Chechk if the second argument is a string
if(~ischar(field))
   error('The second argument ''field'' must be a character string')
end

% cell array 'mn' here is used as comma separated list 
% to pass as an argument to getfield.
switch nargin
   case 5,
      mn = {varargin{1},varargin{2}}; 
   case 4,
      mn = {varargin{1}};
   case 3,
      mn = {1};
   otherwise 
      error('Incorrect number of inputs');
end % switch


NV = length(index);

if(isnumeric(getfield(CELLARRAY{index(1,1)},field)))
   values = zeros(NV,1);
   for I = 1:NV
      values(I) = getfield(CELLARRAY{index(I)},field,mn);
   end 
elseif(ischar(getfield(CELLARRAY{index(1,1)},field)))
   values = cell(NV,1);
   for I = 1:NV
      values{I} = getfield(CELLARRAY{index(I)},field);
   end 
else
   error('The field data must be numeric or character array') 
end



