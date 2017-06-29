function [orbit, fixedpoint] = findsyncorbit(RING,dCT,varargin)
%FINDSYNCORBIT finds closed orbit, synchronous with the RF cavity
% and momentum deviation dP (first 5 components of the phase space vector)
% by numerically solving  for a fixed point
% of the one turn map M calculated with LINEPASS
%
%       (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
%
%    under constraints CT2 - CT1 =  dCT = C(1/Frev - 1/Frev0) and dP2 = dP1 , where
%    Frev0 = Frf0/HarmNumber is the design revolution frequency
%    Frev  = (Frf0 + dFrf)/HarmNumber is the imposed revolution frequency
%
% IMPORTANT!!!  FINDSYNCORBIT imposes a constraint (CT2 - CT1) and
%    dP2 = dP1 but no constraint on the value of dP1, dP2
%    The algorithm assumes time-independent fixed-momentum ring
%    to reduce the dimensionality of the problem.
%    To impose this artificial constraint in FINDSYNCORBIT
%    PassMethod used for any element SHOULD NOT
%    1. change the longitudinal momentum dP (cavities , magnets with radiation)
%    2. have any time dependence (localized impedance, fast kickers etc).
%
%
% FINDSYNCORBIT(RING,dCT) is 5x1 vector - fixed point at the
%		entrance of the 1-st element of the RING (x,px,y,py,dP)
%
% FINDSYNCORBIT(RING,dCT,REFPTS) is 5-by-Length(REFPTS)
%     array of column vectors - fixed points (x,px,y,py,dP)
%     at the entrance of each element indexed in REFPTS array.
%     REFPTS is an array of increasing indexes that  select elements
%     from the range 1 to length(RING)+1.
%     See further explanation of REFPTS in the 'help' for FINDSPOS
%
% FINDSYNCORBIT(RING,dCT,REFPTS,GUESS) - same as above but the search
%     for the fixed point starts at the initial condition GUESS
%     Otherwise the search starts from [0 0 0 0 0 0]'.
%     GUESS must be a 6-by-1 vector;
%
% [ORBIT, FIXEDPOINT] = FINDSYNCORBIT( ... )
%     The optional second return parameter is
%     a 6-by-1 vector of initial conditions
%     on the close orbit at the entrance to the RING.
%
% See also FINDORBIT4, FINDORBIT6.
%
if ~iscell(RING)
    error('First argument must be a cell array');
end

d = 1e-6;       % step size for numerical differentiation
dps = 1e-12;	% convergence threshold
%dps=eps;       % convergence threshold
max_iterations = 20;

if nargin >= 4	% Check if guess argument was supplied
    if isnumeric(varargin{2}) && all(eq(size(varargin{2}),[6,1]))
        Ri=varargin{2};
    else
        error('The last argument GUESS must be a 6-by-1 vector');
    end
else
    Ri = zeros(6,1);
end

D5 = [d*eye(5) zeros(5,1); zeros(1,6)];
%D5 = [0.5*d*eye(5) -0.5*d*eye(5) zeros(5,1); zeros(1,11)];
theta5 = [0 0 0 0  dCT]';

args={};
change=Inf;
itercount = 0;
while (change > dps) && (itercount < max_iterations)
    RMATi= Ri(:,ones(1,6)) + D5;
    %RMATi= Ri(:,ones(1,11)) + D5;
    RMATf = linepass(RING,RMATi,args{:});
    Rf = RMATf(:,end);
    % compute the transverse part of the Jacobian
    J5 =  (RMATf([1:4,6],1:5)-RMATf([1:4,6],6*ones(1,5)))/d;
    %J5 =  (RMATf([1:4,6],1:5)-RMATf([1:4,6],6:10))/d;
    Ri_next = Ri +  [(diag([1 1 1 1 0]) - J5)\(Rf([1:4,6])-[Ri(1:4);0]-theta5);0];
    change = norm(Ri_next - Ri);
    Ri = Ri_next;
    itercount = itercount+1;
    args={'KeepLattice'};
end

if itercount == max_iterations
    warning('Maximum number of iterations reached. Possible non-convergence')
end

if (nargin<3) || (isscalar(varargin{1}) && (varargin{1}==(length(RING)+1)))
    % return only the fixed point at the entrance of RING{1}
    orbit=Ri(1:5,1);
else            % 3-rd input argument - vector of reference points along the Ring
    % is supplied - return orbit
    orb6 = linepass(RING,Ri,varargin{1},'KeepLattice');
    orbit = orb6(1:5,:);
end

if nargout==2
    fixedpoint=Ri;
end
