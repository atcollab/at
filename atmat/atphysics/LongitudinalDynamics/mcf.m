function a = mcf(RING,dp0)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
%
% MCF(RING,DPP) computes the momentum compaction for off-momentum DPP
%
% IMPORTANT!!!
% MCF gives a wrong result with 6-d rings. The RING should be set to 4d.
% See also: ATDISABLE_6D, CHECK_6D


if nargin < 2, dp0=0; end
ddp = 0.000001;
fpdown = findorbit4(RING,dp0-0.5*ddp, 'strict', -1);
fpup = findorbit4(RING,dp0+0.5*ddp, 'strict', -1);
% Build initial condition vector that starts
% on the fixed point

Xdown = [fpdown;dp0-0.5*ddp;0];
Xup = [fpup;dp0+0.5*ddp;0];

% Track X0 for 1 turn
T = ringpass(RING,[Xdown Xup]);
% Calculate alpha
RingLength = findspos(RING,length(RING)+1);
a = (T(6,2)-T(6,1))/(ddp*RingLength);

