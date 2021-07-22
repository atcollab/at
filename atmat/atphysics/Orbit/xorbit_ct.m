function Ri = xorbit_ct(ring,dct,varargin)
%XORBIT_CT	Private function used by findsyncorbit
%ORBIT=XORBIT_CT(RING,DCT)	Closed orbit for fixed path lengthening.
[XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation
[dps,varargs]=getoption(varargs,'OrbConvergence');	% Convergence threshold
[max_iterations,varargs]=getoption(varargs,'OrbMaxIter');	% Max. iterations
[Ri,varargs]=getoption(varargs,'guess',zeros(6,1));
[Ri,varargs]=getargs(varargs,Ri,'check',@(arg) isnumeric(arg) && isequal(size(arg),[6,1])); %#ok<ASGLU>

if ~isfinite(dct), dct = 0.0; end
scaling=XYStep*[1 1 1 1 1];
D5 = [diag(scaling) zeros(5,1); zeros(1,6)];
theta5 = [0 0 0 0  dct]';

change=Inf;
itercount = 0;
args={};
while (change > dps) && (itercount < max_iterations)
    RMATi= Ri(:,ones(1,6)) + D5;
    RMATf = linepass(ring,RMATi,args{:});
    Rf = RMATf(:,end);
    % compute the transverse part of the Jacobian
    J5 =  (RMATf([1:4,6],1:5)-RMATf([1:4,6],6*ones(1,5)))./scaling;
    Ri_next = Ri +  [(diag([1 1 1 1 0]) - J5)\(Rf([1:4,6])-[Ri(1:4);0]-theta5);0];
    change = norm(Ri_next - Ri);
    Ri = Ri_next;
    itercount = itercount+1;
    args={'KeepLattice'};
end

if itercount == max_iterations
    warning('Maximum number of iterations reached. Possible non-convergence')
end
end
