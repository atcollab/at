function spos = findspos(LINE,REFPTS)
%FINDSPOS returns longitudinal positions of accelerator lattice elements.
%  Return value is a vector of positions S at the entrance of each element
%  specified by its number in REFPTS
%
% Note: REFPTS is an array of increasing indexes that  
%   select elements from range 1 to length(LATTICE)+1. 
%   REFPTS is allowed to go 1 point beyond the 
%   number of elements. In this case the last point is 
%   the EXIT of the last element. If LATTICE is a RING
%   it is also the entrance of the first element after 1 turn.
%    
%   Note:
%   1. Use findspos(RING,1:length(RING)) for to find 
%      longitudinal position of all elements 
%   2. Use findspos(LINE,length(LINE)+1) to find the 
%      total physical length
%   3. If line is a closed ring, exit of the last element 
%      is also the entrance to the first. 
    
NE=length(LINE);
L=zeros(1,NE+1);

for k=2:NE+1
   L(k)=L(k-1)+LINE{k-1}.Length;
end;
spos=L(REFPTS);