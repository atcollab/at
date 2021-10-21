function Ri = xorbit_6(RING,varargin)
%XORBIT_6	Private function used by findorbit6
%ORBIT=XORBIT_6(RING)	Closed orbit for 6D motion
[XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation	%1.e-8
[DPStep,varargs]=getoption(varargs,'DPStep');	% Step size for numerical differentiation	%1.e-6
[dps,varargs]=getoption(varargs,'OrbConvergence');	% Convergence threshold                 %1.e-12
[max_iterations,varargs]=getoption(varargs,'OrbMaxIter');	% Max. iterations               %20
[method,varargs]=getoption(varargs,'method','tracking');    % Method for Eloss computation
[Ri,varargs]=getoption(varargs,'guess',[]);
[Ri,varargs]=getargs(varargs,Ri,'check',@(arg) isnumeric(arg) && isequal(size(arg),[6,1])); %#ok<ASGLU>

cavities=RING(atgetcells(RING,'Frequency'));
if isempty(cavities)
    error('AT:NoFrequency', 'The lattice has no Cavity element')
end
Frf=min(atgetfieldvalues(cavities,'Frequency'));

L0 = findspos(RING,length(RING)+1); % design length [m]
C0 = PhysConstant.speed_of_light_in_vacuum.value; % speed of light [m/s]
T0 = L0/C0;                         % Revolution period [s]
HarmNumber = round(Frf*L0/C0);

if isempty(Ri)
    Vcell=sum(atgetfieldvalues(cavities,'Voltage'));
    U0cell=atgetU0(RING,'periods',1,'method',method);
    if U0cell > Vcell
        error('AT:MissingVoltage','Missing RF voltage: unstable ring');
    end
    Ri = zeros(6,1);
    Ri(6) = -L0/(2*pi*HarmNumber) * asin(U0cell/Vcell);
end

theta = [0 0 0 0 0 C0*(HarmNumber/Frf - T0)]';

scaling=XYStep*[1 1 1 1 0 0] + DPStep*[0 0 0 0 1 1];
D = [diag(scaling) zeros(6,1)];

change=Inf;
itercount = 0;
args={};
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
end

