function Ri = xorbit_dp(ring,dp,varargin)
%XORBIT_DP  Private function used by findorbit4
%ORBIT=XORBIT_DP(RING,DP)   Closed orbit for fixed off-momentum.
[XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation
[dps,varargs]=getoption(varargs,'OrbConvergence');	% Convergence threshold
[max_iterations,varargs]=getoption(varargs,'OrbMaxIter');	% Max. iterations
[Ri,varargs]=getoption(varargs,'guess',zeros(6,1));
[Ri,varargs]=getargs(varargs,Ri,'check',@(arg) isnumeric(arg) && isequal(size(arg),[6,1])); %#ok<ASGLU>

% Set the momentum component of Ri to the specified dP
if ~isfinite(dp), dp = 0.0; end
Ri(5) = dp;
scaling=XYStep*[1 1 1 1];
D = [diag(scaling) zeros(4,1);zeros(2,5)];

change=Inf;
itercount = 0;
args={};
while (change > dps) && (itercount < max_iterations)
    RMATi = Ri(:,ones(1,5)) + D;
    RMATf = linepass(ring,RMATi,args{:});
    Rf = RMATf(:,end);
    % compute the transverse part of the Jacobian
    J4 = (RMATf(1:4,1:4)-RMATf(1:4,5))./scaling;
    Ri_next = Ri +  [(eye(4) - J4)\(Rf(1:4)-Ri(1:4)); 0; 0];
    change = norm(Ri_next - Ri);
    Ri = Ri_next;
    itercount = itercount+1;
    args={'KeepLattice'};
end

if itercount == max_iterations
    warning('Maximum number of iterations reached. Possible non-convergence')
end
end
