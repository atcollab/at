function opts = initoptions()
%INITOPTIONS	Set initial values for AT constants
%
%OPTIONS=INITOPTIONS()
%OPTIONS: Scalar structure with predefined default values for AT constants

opts.XYStep = 1.0E-8;  % 6.055454452393343e-006
opts.DPStep = 1.0E-6;
opts.OrbConvergence = 1.0E-12;
opts.OrbMaxIter = 20;
end

