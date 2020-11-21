function [orb4,fixedpoint] = findorbit4(RING,dP,varargin)
%FINDORBIT4 finds closed orbit in the 4-d transverse phase
% space by numerically solving  for a fixed point of the one turn
% map M calculated with LINEPASS
%
%         (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
%
%    under the CONSTANT MOMENTUM constraint, dP2 = dP1 = dP and
%    there is NO constraint on the 6-th coordinate CT
%
% IMPORTANT!!! FINDORBIT4 imposes a constraint on dP and relaxes
%    the constraint on the revolution frequency. A physical storage
%    ring does exactly the opposite: the momentum deviation of a
%    particle on the closed orbit settles at the value
%    such that the revolution is synchronous with the RF cavity
%
%                 HarmNumber*Frev = Frf
%
%    To impose this artificial constraint in FINDORBIT4
%    PassMethod used for any elemen SHOULD NOT
%    1. change the longitudinal momentum dP (cavities , magnets with radiation)
%    2. have any time dependence (localized impedance, fast kickers etc)
%
% FINDORBIT4(RING,dP) is 4x1 vector - fixed point at the
%    entrance of the 1-st element of the RING (x,px,y,py)
%
% FINDORBIT4(RING,dP,REFPTS) is 4-by-Length(REFPTS)
%     array of column vectors - fixed points (x,px,y,py)
%     at the entrance of each element indexed  REFPTS array.
%     REFPTS is an array of increasing indexes that  select elements
%     from the range 1 to length(RING)+1.
%     See further explanation of REFPTS in the 'help' for FINDSPOS
%
% FINDORBIT4(RING,dP,REFPTS,GUESS) - same as above but the search
%     for the fixed point starts at the initial condition GUESS
%     Otherwise the search starts from [0 0 0 0 0 0]'.
%     GUESS must be a 6-by-1 vector;
%
% [ORBIT, FIXEDPOINT] = FINDORBIT4( ... )
%     The optional second return parameter is
%     a 6-by-1 vector of initial conditions
%     on the closed orbit at the entrance to the RING.
%
% See also FINDSYNCORBIT, FINDORBIT6.

if ~iscell(RING)
    error('First argument must be a cell array');
end
[XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation
[dps,varargs]=getoption(varargs,'OrbConvergence');	% Convergence threshold
[max_iterations,varargs]=getoption(varargs,'OrbMaxIter');	% Max. iterations

if length(varargs) >= 1 && ~isequal(varargs{1},length(RING)+1)
    refpts=varargs{1};
    if islogical(refpts)
        refpts=find(refpts);
    end
else
    refpts=[];
end
if length(varargs) >= 2	% Check if guess argument was supplied
    if isnumeric(varargs{2}) && isequal(size(varargs{2}),[6,1])
        Ri=varargs{2};
    else
        error('The last argument GUESS must be a 6-by-1 vector');
    end
else
    Ri = zeros(6,1);
end

% Set the momentum component of Ri to the specified dP
Ri(5) = dP;
D = [XYStep*eye(4) zeros(4,1);zeros(2,5)];
%D = [0.5*d*eye(4) -0.5*d*eye(4) zeros(4,1);zeros(2,9)];

args={};
change=Inf;
itercount = 0;
while (change > dps) && (itercount < max_iterations)
    RMATi = Ri(:,ones(1,5)) + D;
    %RMATi = Ri(:,ones(1,9)) + D;
    RMATf = linepass(RING,RMATi,args{:});
    Rf = RMATf(:,end);
    % compute the transverse part of the Jacobian
    J4 = (RMATf(1:4,1:4)-RMATf(1:4,5*ones(1,4)))/XYStep;
    %J4 = (RMATf(1:4,1:4)-RMATf(1:4,5:8))/d;
    Ri_next = Ri +  [(eye(4) - J4)\(Rf(1:4)-Ri(1:4)); 0; 0];
    change = norm(Ri_next - Ri);
    Ri = Ri_next;
    itercount = itercount+1;
    args={'KeepLattice'};
end

if itercount == max_iterations
    warning('Maximum number of iterations reached. Possible non-convergence')
end

if isempty(refpts)
    % return only the fixed point at the entrance of RING{1}
    orb4=Ri(1:4,1);
else	% 3-rd input argument - vector of reference points along the RING
    % is supplied - return orbit
    orb6 = linepass(RING,Ri,refpts,'KeepLattice');
    orb4 = orb6(1:4,:);
end

if nargout >= 2
    fixedpoint=Ri;
end
