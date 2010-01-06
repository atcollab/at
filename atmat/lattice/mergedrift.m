function mergedrift(SPLITPOS) 
%MERGEDRIFT removes a lattice element and merges the two adjacent drift spaces
%
% MERGEDRIFT (SPLITPOS) removes an element located at SPLITPOS from the global lattice THERING
% surrounded by two DRIFT spaces. The resulting drift has Length L0 = L1 + LSPLIT + L2;  
% Number of elements in THERING is thus reduced by 2
%
% See also: SPLITDRIFT

global THERING
L0 = THERING{SPLITPOS-1}.Length  + THERING{SPLITPOS+1}.Length + THERING{SPLITPOS}.Length;


% make new (empty) THERING 
NEWN = length(THERING)-2;
R = cell(1,NEWN);

N0 = length(THERING); % number of elements in THERING before split

for i = 1:(SPLITPOS-1)
	R{i} = THERING{i};
end

R{SPLITPOS-1}.Length = L0;

for i = SPLITPOS:NEWN
	R{i} = THERING{i+2};
end

THERING = R;
