function mergedrift(ELEMPOS, DIST) 
%MVELEM(ELEMPOS, DIST) moves an element  located at ELEMPOS in THERING
% surrounded by two DRIFT spaces 
% 
%  0   < DIST  < LD move downstream
% -LU  < DIST  < 0  move upstream
%  where LU and LD - lenths of 
%  upstream and downstrem drift drifts BEFORE!!! the move 
%
% Number of elements in THERING and total length remain the same
%
% See also: SPLITDRIFT, MERGEDRIFT

global THERING

L0 = THERING{ELEMPOS-1}.Length + THERING{ELEMPOS}.Length + THERING{ELEMPOS+1}.Length;

if DIST > THERING{ELEMPOS+1}.Length 
	error('Cannot move downstream more than the length of downstream drift');
elseif -DIST > THERING{ELEMPOS-1}.Length 
	error('Cannot move upstream more than the length of upstream drift');
else
	THERING{ELEMPOS+1}.Length = THERING{ELEMPOS+1}.Length - DIST;
	THERING{ELEMPOS-1}.Length = THERING{ELEMPOS-1}.Length + DIST;
end

		
	 
