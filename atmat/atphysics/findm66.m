function [M66, varargout] = findm66(RING, varargin);
%FINDM66 numerically finds the 6x6 transfer matrix of an accelerator lattice
%  by differentiation of LINEPASS near the closed orbit
%  FINDM66 uses fFINDORBIT6 to search for the closed orbit in 6-d
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
%    When REFPTS is a vector FINDM44 is a 6-by-6-by-length(REFPTS) array
%    so that the set of indexes (:,:,i) selects the 6-by-6 
%    matrix at the i-th reference point
%    
%    Note: 
%    When REFPTS= [ 1 2 .. ] the fist point is the entrance of the first element
%    and T(:,:,1) - identity matrix
%    When REFPTS= [  .. length(RING)+1] the last point is the exit of the last element
%    and the entrance of the first element after 1 turn: T(:,:, ) = M
%
% [M66, T, orbit] = FINDM66(RING, REFPTS) in addition returns the closed orbit
%    found in the process of lenearization
%
% See also findm44

FULL = 0;
switch nargin 
    
    case 2 
      REFPTS = varargin{1};
   case 1
   	REFPTS = 1;
	otherwise
     error('Too many input arguments')
end


NE = length(RING);
NR = length(REFPTS);

% See if step size for numerical differentiation
% is set globally. Otherwise use 1e-7
global NUMDIFPARAMS
% Transverse
if isfield(NUMDIFPARAMS,'XYStep')
    dt = NUMDIFPARAMS.XYStep';
else
    dt =  1e-8;
end
% Longitudinal
if isfield(NUMDIFPARAMS,'DPStep')
    dl = NUMDIFPARAMS.DPStep';
else
    dl =  1e-8;
end


% Calculate closed orbit in 6 dimensions (MUST have cavity in the ring)
reforbit = findorbit6(RING,REFPTS);
% Build a diagonal matrix of initial conditions
D6 = [dt*eye(4),zeros(4,2);zeros(2,4), dl*eye(2)];
% Add to the orbit_in
RIN = reforbit(:,1)*ones(1,12) + [D6, -D6];




if nargout <= 1 % Whole ring , NO REFPTS
    % Propagate through the ring
    ROUT = ringpass(RING,RIN);
    % Calculate numerical derivative
    M66 = [(ROUT(:,1:4)-ROUT(:,7:10))./(2*dt), (ROUT(:,5:6)-ROUT(:,11:12))./(2*dl)];
   return
else					
    % Calculate matrixes at all REFPTS. Use linepass
    % Need to include the exit of the RING to REFPTS array
    if(REFPTS(NR)~=NE+1)
        REFPTS = [REFPTS NE+1];
        NR1 = NR+1;
    else
        NR1 = NR;
    end
	TMAT = linepass(RING,RIN,REFPTS);
	% Reshape, so that the otcome at each REFPTS is a separate page in a 3-dim array
    TMAT3 = reshape(TMAT,6,12,NR1);
    varargout{1} = [(TMAT3(:,1:4,1:NR)-TMAT3(:,7:10,1:NR))/(2*dt), (TMAT3(:,5:6,1:NR)-TMAT3(:,11:12,1:NR))/(2*dl)];
    M66 = [(TMAT3(:,1:4,NR1)-TMAT3(:,7:10,NR1))/(2*dt), (TMAT3(:,5:6,NR1)-TMAT3(:,11:12,NR1))/(2*dl)];
end

% Return closed orbit if requested
if nargout == 3
   varargout{2} = reforbit;
end
