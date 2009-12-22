function index = findtags(CELLARRAY, MATCHSTR)
%FINDTAGS looks for string matches in 'Tag' field of AT lattice elements
%   
% INDEX = FINDTAGS(CELLARRAY, MATCHSTR) 
%   returns indexes of elements that have a field 'Tag'
%   whose value is a string exactly matching MATCHSTR
%   or a cell array of strings with one element matching MATCHSTR
%
% See also FINDCELLS, SETTAGS, 

% Check if the first argument is the cell arrray of tstructures 
if(~iscell(CELLARRAY) | ~isstruct(CELLARRAY{1}) | isempty(CELLARRAY))
   error('The first argument must be a non-empty cell array of structures') 
end
% Chechk if the second argument is a string
if(~ischar(MATCHSTR))
      error('The second argument must be a character string')
end




NE = length(CELLARRAY);
matchesfound = 0;
index = findcells(CELLARRAY,'Tag');


index1 = index;
matchesfound = 0;
for I = index
    if any(strcmp(CELLARRAY{I}.Tag,MATCHSTR));
        matchesfound = matchesfound+1;
        % since 'matchesfound' counter is <= loop number,
        % it is save to modify elements of 'index' inside the loop
        index(matchesfound) = I;

    end
end

index =  index(1:matchesfound);

        
