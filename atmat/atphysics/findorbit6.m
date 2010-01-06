function orbit = findorbit6(RING,varargin);
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
% See also FINDORBIT, FINDORBIT4, FINDSYNCORBIT.

if ~iscell(RING)
   error('First argument must be a cell array');
end 


L0 =   findspos(RING,length(RING)+1); % design length [m]
C0 =   299792458; % speed of light [m/s]
T0 = L0/C0;

CavityLocation = findcells(RING,'Frequency');
Frf = RING{CavityLocation(1)}.Frequency;

if ~isfield(RING{CavityLocation(1)},'HarmNumber')
    error('FINDORBIT6 needs the HarmNumber field to be defined in RF cavity elements');
end
HarmNumber = RING{CavityLocation(1)}.HarmNumber;
theta = [0 0 0 0 0 C0*(HarmNumber/Frf - T0)]';

 
d = 1e-6;	% step size for numerical differentiation
max_iterations = 20;

if nargin==3
    if isnumeric(varargin{2}) & all(eq(size(varargin{2}),[6,1]))
       Ri=varargin{2};
   else
       error('The last argument GUESS must be a 6-by-1 vector');
   end
else
    Ri = zeros(6,1);
end
D = d*eye(6);  

RMATi= [Ri Ri Ri Ri Ri Ri Ri] + [D, zeros(6,1)];
RMATf = linepass(RING,RMATi);
J6 = (RMATf(:,1:6)-RMATf(:,7)*ones(1,6))/d;
Rf = RMATf(:,7);
% Replace matrix inversion with \
% B = inv(eye(6)-J6);
% Ri_next = Ri + B*(Rf-Ri-theta);
Ri_next = Ri + (eye(6)-J6)\(Rf-Ri-theta);
change = norm(Ri_next - Ri);
Ri = Ri_next;

itercount = 1;


while (change>eps) & (itercount < max_iterations)
   RMATi= [Ri Ri Ri Ri Ri Ri Ri] + [D, zeros(6,1)];
   RMATf = linepass(RING,RMATi,'reuse');
   J6 = (RMATf(:,1:6)-RMATf(:,7)*ones(1,6))/d;
   Rf = RMATf(:,7);
% Replace matrix inversion with \
%    B = inv(eye(6)-J6);
%    Ri_next = Ri + B*(Rf-Ri-theta);
    Ri_next = Ri + (eye(6)-J6)\(Rf-Ri-theta);
   change = norm(Ri_next - Ri);
   Ri = Ri_next;
   itercount = itercount+1;
end;

if itercount == max_iterations
    warning('Maximum number of itereations reached. Possible non-convergence')
end
if(nargin==1)|(varargin{1}==(length(RING)+1))
   % return only the fixed point at the entrance of RING{1}
   orbit=Ri;
else % 2-nd input argument - vector of reference points alog the Ring
                  % is supplied - return orbit            
	orbit = linepass(RING,Ri,varargin{1},'reuse'); 
end

if nargout==2
    varargout{1}=Ri;
end