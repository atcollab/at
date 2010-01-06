function a = mcf(RING)
%MCF momentum compaction factor
% MCF(RING) calculates momentum compaction factor of RING
dP = 0.000001;
fp0 = findorbit4(RING,0);
fp = findorbit4(RING,dP);
% Build initial condition vector that starts
% on the fixed point
X0dP = fp;
X0dP(5) = dP;
X0dP(6) = 0;

X0 = [fp0;0;0];

% Track X0 for 1 turn
T = ringpass(RING,[X0 X0dP]);
% Calculate alpha
RingLength = findspos(RING,length(RING)+1);
a = (T(6,2)-T(6,1))/(dP*RingLength);

