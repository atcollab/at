function [M66, varargout] = findm66(LATTICE, varargin)
%FINDM66 numerically finds the 6x6 transfer matrix of an accelerator lattice
%  by differentiation of LINEPASS near the closed orbit
%  FINDM66 uses FINDORBIT6 to search for the closed orbit in 6-d
%  In order for this to work the ring MUST have a CAVITY element
%
% M66 = FINDM66(RING)finds full one-turn 6-by-6
%    matrix at the entrance of the first element
%
% [M66,T] = FINDM66(RING,REFPTS) in addition to M finds
%    6-by-6 transfer matrixes  between entrances of
%    the first element and each element indexed by REFPTS.
%    T is 6-by-6-by-length(REFPTS) 3 dimentional array.
%
%    REFPTS is an array of increasing indexes that  select elements
%    from the range 1 to length(RING)+1.
%    See further explanation of REFPTS in the 'help' for FINDSPOS
%
%    Note:
%    When REFPTS= [ 1 2 .. ] the fist point is the entrance of the first element
%    and T(:,:,1) - identity matrix
%    When REFPTS= [  .. length(RING)+1] the last point is the exit of the last element
%    and the entrance of the first element after 1 turn: T(:,:, ) = M
%
% [M66,T] = FINDM66(RING,REFPTS,ORBITIN) - Does not search for closed orbit.
% [...]   = FINDM66(...,'orbit',ORBITIN)
%   This syntax is useful to avoid recomputing the closed orbit if it is
%   already known.
%
% [M66,T,orbit] = FINDM66(RING, REFPTS) in addition returns the closed orbit
%    found in the process of lenearization
%
% See also FINDM44, FINDORBIT6

if ~iscell(LATTICE)
    error('First argument must be a cell array');
end
NE = length(LATTICE);
[XYStep,varargs]=getoption(varargin,'XYStep');
[R0,varargs]=getoption(varargs,'orbit',missing);

if length(varargs) >= 2	% FINDM66(RING,REFPTS,ORBITIN)
    R0 = varargs{2};
end

if  ismissing(R0)
    R0 = findorbit6(LATTICE);
end

if length(varargs) >= 1	% FINDM66(RING,REFPTS)
    if islogical(varargs{1})
        REFPTS=varargs{1};
        REFPTS(end+1:NE+1)=false;
    elseif isnumeric(varargs{1})
        REFPTS=setelems(false(1,NE+1),varargs{1});
    else
        error('REFPTS must be numeric or logical');
    end
else
    REFPTS=false(1,NE+1);
end
refs=setelems(REFPTS,NE+1);
reqs=REFPTS(refs);

% Build a diagonal matrix of initial conditions
%scaling=2*XYStep*[1 0.1 1 0.1 1 1];
scaling=2*XYStep*[1 1 1 1 1 1];
D6 = 0.5*diag(scaling);
% Add to the orbit_in. First 12 columns for derivative
% 13-th column is for closed orbit
RIN = R0 + [D6 -D6 zeros(6,1)];
ROUT = linepass(LATTICE,RIN,refs);
TMAT3 = reshape(ROUT,6,13,[]);
M66 = (TMAT3(:,1:6,end)-TMAT3(:,7:12,end))./scaling;

if nargout >= 2 % Calculate matrices at all REFPTS.
    varargout{1} = (TMAT3(:,1:6,reqs)-TMAT3(:,7:12,reqs))./scaling;
    % Return closed orbit if requested
    if nargout >= 3
        varargout{2}=squeeze(TMAT3(:,13,reqs));
    end
end

    function mask=setelems(mask,idx)
        mask(idx)=true;
    end

end
