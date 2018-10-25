function [ t ] = twiss_transfer(RING,DP,refpts)
twiss = twissring(RING,DP,refpts);
Orbit = [twiss.ClosedOrbit];
M44 = reshape([twiss.M44], 4, 4, []);
Alpha = reshape([twiss.alpha], 2, []);
Beta = reshape([twiss.beta], 2, []);
Mu = reshape([twiss.mu], 2, []);
t = struct('O', Orbit, 'M', M44, 'A', Alpha, 'B', Beta, 'U', Mu);
end