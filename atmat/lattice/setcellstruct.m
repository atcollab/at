function CELLARRAY = setcellstruct(CELLARRAY,field,index,value,varargin)
%SETCELLSTRUCT sets the field values of MATLAB cell array of structures
% 
%   Note that the calling syntax must be in the form of assignment:
%   CELLARRAY = SETCELLSTRUCT(CELLARRAY,...)
%   MATLAB does not modify variables that only appear on the right 
%   hand side as arguments. 
%
% Numeric data
% ---------------------------------------------------------
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M,N)
%	 Sets (M,N) element equal to VALUE when the field data is 
%   a mtrix. The assigned VALUE may be 
%   1. Numeric array of the same length as INDEX array
%   2. Single numeric value to be written to CELLARRAY elements
%
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M) can be used 
%   for one dimensional vectors
%
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is the same as 
%   SETCELLSTRUCT(CELLARRAY,'field',index,1,1) if the field data
%   is a scalar
%
% Character array 
% --------------------------------------------------------------------
% CELLARRAY SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is a MATLAB 
%   cell array of strings when specified fields contain strings.
%   The assignment VALUE may be
%   1. Character array with the number of rows matching the number of 
%      elements in INDEX array
%   2. Character string
%
% See also GETCELLSTRUCT FINDCELLS
global THERING 
if (nargout<1)
   error('Must have output arguments. Use assignment syntax')
end
    
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
   case 6,
      mn = {varargin{1},varargin{2}}; 
   case 5,
      mn = {varargin{1}};
   case 4,
      mn = {1};
   otherwise 
      error('Incorrect number of inputs');
end % switch

vtype = class(value);
switch vtype
   case 'double',
      if(length(index)==length(value))
         for I=1:length(index)
            CELLARRAY{index(I)} = setfield(CELLARRAY{index(I)},field,mn,value(I));
         end
      elseif(length(value)==1)
         for I=1:length(index)
            CELLARRAY{index(I)} = setfield(CELLARRAY{index(I)},field,mn,value);
         end
      else
         error('Number of elements in VALUE must be 1 or match the length of INDEX');
      end
   case 'char',
      if(length(index)==size(value,1))
         for I=1:length(index)
            CELLARRAY{index(I)} = setfield(CELLARRAY{index(I)},field,value(I,:));
         end
      elseif(size(value,1)==1)
         for I=1:length(index)
            CELLARRAY{index(I)} = setfield(CELLARRAY{index(I)},field,value);
         end
      else
         error('Number of rows in character array VALUE must be 1 or match the length of INDEX');
      end
   otherwise
      error('The field data must be numeric or character array');
 end
 


   
