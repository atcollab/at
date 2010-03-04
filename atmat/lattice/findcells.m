function index = findcells(CELLARRAY, field, varargin)
%FINDCELLS performs a search on MATLAB cell arrays of structures
%   
% INDEX = FINDCELLS(CELLARRAY, 'field') 
%   returns indexes of elements that have a field named 'field'   
%
% INDEX = FINDCELLS(CELLARRAY, 'field', VALUE) 
%   returns indexes of elements whose field 'field'
%   is equal to VALUE1, VALUE2, ... or VALUEN. Where VALUE can either be
%   character strings or a number. If its a character string REGULAR
%   expressions can be used.
%
% Example:
%   findcells(THERING,'Length',0, 0.2);  % will match elements of
%                                          lengths 0 and 0.2
%   findcells(THERING,'FamName','SFA','SDA');
%
% See also GETCELLSTRUCT, SETCELLSTRUCT, REGEXPI

% Check if the first argument is the cell arrray of tstructures 
if(~iscell(CELLARRAY) | ~isstruct(CELLARRAY{1}) | isempty(CELLARRAY))
   error('The first argument must be a non-empty cell array of structures') 
end
% Chechk if the second argument is a string
if(~ischar(field))
      error('The second argument must be a character string')
end
% if(nargin > 3)
%      error('Incorrect number of inputs')
% end



NE = length(CELLARRAY);
matchesfound = 0;
index = zeros(1,NE);
for I = 1:NE
   if(isfield(CELLARRAY{I},field))
      matchesfound = matchesfound+1;
      index(matchesfound) = I; 
   end
end

index =  index(1:matchesfound);
if(nargin > 2)
    index1 = index;
    matchesfound = 0;
    for I = index
        for j=1:length(varargin)
            if ~isnumeric(varargin{j})
                %if isequal(getfield(CELLARRAY{I},field),varargin{j})
                if regexpi(getfield(CELLARRAY{I},field),['^' varargin{j} '$'])
                    matchesfound = matchesfound+1;
                    % since 'matchesfound' counter is <= loop number,
                    % it is save to modify elements of 'index' inside the loop
                    index(matchesfound) = I;
                    
                end
            else
                %if you are comparing a numerical value
                if getfield(CELLARRAY{I},field) == varargin{j}
                    matchesfound = matchesfound+1;
                    % since 'matchesfound' counter is <= loop number,
                    % it is save to modify elements of 'index' inside the loop
                    index(matchesfound) = I;
                end
            end
        end
    end
    
    index =  index(1:matchesfound);
end


        
