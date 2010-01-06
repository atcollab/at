function splitdrift(DRIFTPOS, SPLIT, varargin) 
%SPLITDRIFT inserts an element into a drift space
%
% SPLITDRIFT(DRIFTPOS, SPLIT) inserts a marker (zero-length) element
%   at distance SPLIT ( 0 < SPLIT < 1) into a drift space 
%   located at DRIFTPOS in THERING
% 
% SPLITDRIFT(DRIFTPOS, SPLIT, ELEMSTRUCCTURE) inserts a marker (zero-length) element
%   at distance SPLIT ( 0 < SPLIT < 1) into a drift space 
%   located at DRIFTPOS in THERING
% 
% Number of elements in the RING is thus increased by 2
% SPLIT (controls the position of the split 
% L1 = L0*SPLIT
% L2 = L0(1-SPLIT)
%  where L0 is the length of the original DRIFT
%   
% See also: MERGEDRIFT
 
global THERING

N0 = length(THERING); % number of elements in THERING before split

L0 = THERING{DRIFTPOS}.Length;
if (SPLIT < 0 | SPLIT >1)
	error('Second argument must be (0..1)');
else	
	if nargin == 2
        for i = reverse(DRIFTPOS:N0)
	        THERING{i+2} = THERING{i};
        end
        THERING{DRIFTPOS}.Length = L0*SPLIT;
	    THERING{DRIFTPOS+2}.Length = L0*(1-SPLIT);
	    THERING{DRIFTPOS+1} = struct('FamName','TEMPSPLIT', 'Length', 0, 'PassMethod','IdentityPass');
    elseif isstruct(varargin{1})
         % Check if a new element will fit
        if L0*(1-SPLIT) >= varargin{1}.Length
            for i = reverse(DRIFTPOS:N0)
	            THERING{i+2} = THERING{i};
            end
            THERING{DRIFTPOS}.Length = L0*SPLIT;
	        THERING{DRIFTPOS+2}.Length = L0*(1-SPLIT)-varargin{1}.Length;
            THERING{DRIFTPOS+1} = varargin{1};
        else
            error('The inserted element is too long');
        end
    end

end
