function spos = findspos(line,refpts)
%FINDSPOS returns longitudinal positions of accelerator lattice elements.
%  Return value is a row vector of positions S at the entrance of each
%  element specified REFPTS (index list or logical mask)
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

if isvector(line)
    L=[0;cumsum(cellfun(@(el) el.Length,line(:),'ErrorHandler',@(w,el) 0))];
else
    L=[0;cumsum(cellfun(@(el) el.Length,line(:,1),'ErrorHandler',@(w,el) 0))];
end
spos=L(refpts)';
end
