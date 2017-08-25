function a = mcf(RING,dp0)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
%
% MCF(RING,DPP) computes the momentum compaction for off-momentum DPP

if nargin < 2, dp0=0; end
ddp = 0.000001;
fp0 = findorbit4(RING,dp0);
fp = findorbit4(RING,dp0+ddp);
% Build initial condition vector that starts
% on the fixed point
X0dP = fp;
X0dP(5) = ddp;
X0dP(6) = 0;

X0 = [fp0;0;0];

% Track X0 for 1 turn
T = ringpass(RING,[X0 X0dP]);
% Calculate alpha
RingLength = findspos(RING,length(RING)+1);
a = (T(6,2)-T(6,1))/(ddp*RingLength);

