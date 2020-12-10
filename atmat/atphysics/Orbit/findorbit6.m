function [orb6,fixedpoint] = findorbit6(RING,varargin)
%FINDORBIT6 finds closed orbit in the full 6-d phase space
% by numerically solving  for a fixed point of the one turn
% map M calculated with RINGPASS
%
% (X, PX, Y, PY, DP, CT2 ) = M (X, PX, Y, PY, DP, CT1)
%
% with constraint % CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
%
% IMPORTANT!!! FINDORBIT6 is a realistic simulation
% 1. The Frf frequency in the RF cavities (may be different from Frf0)
%    imposes the synchronous condition
%    CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
% 2. The algorithm numerically calculates
%    6-by-6 Jacobian matrix J6. In order for (J-E) matrix
%    to be non-singular it is NECESSARY to use a realistic
%    PassMethod for cavities with non-zero momentum kick
%    (such as ThinCavityPass).
% 3. FINDORBIT6 can find orbits with radiation.
%    In order for the solution to exist the cavity must supply
%    adequate energy compensation.
%    In the simplest case of a single cavity, it must have
%    'Voltage' field set so that Voltage > Erad - energy loss per turn
% 4. FINDORBIT6 starts the search from [ 0 0 0 0 0 0 ]', unless
%    the third argument is specified: FINDORBIT6(RING,REFPTS,GUESS)
%    There exist a family of solutions that correspond to different RF buckets
%    They differ in the 6-th coordinate by C*Nb/Frf. Nb = 1 .. HarmNum-1
% 5. The value of the 6-th coordinate found at the cavity gives
%    the equilibrium RF phase. If there is no radiation the phase is 0;
%
% FINDORBIT6(RING) is 6x1 vector - fixed point at the
%		entrance of the 1-st element of the RING (x,px,y,py,dp,ct)
%
% FINDORBIT6(RING,REFPTS) is 6-by-Length(REFPTS)
%     array of column vectors - fixed points (x,px,y,py,dp,ct)
%     at the entrance of each element indexed  REFPTS array.
%     REFPTS is an array of increasing indexes that  select elements
%     from the range 1 to length(RING)+1.
%     See further explanation of REFPTS in the 'help' for FINDSPO
%
% FINDORBIT6(RING,REFPTS,GUESS) - same as above but the search
%     for the fixed point starts at the initial condition GUESS
%     GUESS must be a 6-by-1 vector;
%
% [ORBIT, FIXEDPOINT] = FINDORBIT6( ... )
%     The optional second return parameter is
%     a 6-by-1 vector of initial conditions
%     on the close orbit at the entrance to the RING.
%
% See also FINDORBIT4, FINDSYNCORBIT.

if ~iscell(RING)
    error('First argument must be a cell array');
end
[XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation	%1.e-6
[dps,varargs]=getoption(varargs,'OrbConvergence');	% Convergence threshold                 %1.e-12
[max_iterations,varargs]=getoption(varargs,'OrbMaxIter');	% Max. iterations               %20

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

L0 = findspos(RING,length(RING)+1); % design length [m]
C0 = PhysConstant.speed_of_light_in_vacuum.value; % speed of light [m/s]

T0 = L0/C0;

CavityIndex = find(atgetcells(RING,'Frequency'),1);

if isempty(CavityIndex)
    error('findorbit6: The lattice does not have Cavity element')
end

cavity1=RING{CavityIndex};

Frf = cavity1.Frequency;
HarmNumber = cavity1.HarmNumber;
theta = [0 0 0 0 0 C0*(HarmNumber/Frf - T0)]';

scaling=XYStep*[1 1 1 1 1 1];
D = [diag(scaling) zeros(6,1)];

args={};
change=Inf;
itercount = 0;
while (change > dps) && (itercount < max_iterations)
    RMATi= Ri(:,ones(1,7)) + D;
    RMATf = linepass(RING,RMATi,args{:});
    % compute the transverse part of the Jacobian 
    J6 = (RMATf(:,1:6)-RMATf(:,7))./scaling;
    Rf = RMATf(:,end);
    Ri_next = Ri + (eye(6)-J6)\(Rf-Ri-theta);
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
    orb6=Ri;
else	% 2-nd input argument - vector of reference points along the Ring
        % is supplied - return orbit
    orb6 = linepass(RING,Ri,refpts,'KeepLattice');
end

if nargout >= 2
    fixedpoint=Ri;
end
