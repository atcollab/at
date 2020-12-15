function opts = initoptions()
%INITOPTIONS	Set initial values for AT constants
%
%OPTIONS=INITOPTIONS()
%
%OPTIONS: Scalar structure with predefined default values for AT constants:
%  XYStep = 3.0E-8;    Transverse step (for findorbitx, findmxx,...)
%  DPStep = 3.0E-6;    Momentum step (for dispersion and chromaticity)
%  OrbConvergence = 1.0E-12;   Convergence of findorbitx
%  OrbMaxIter = 20;            Max number of iterations for findorbitx
%
%see also GETOPTION, SETOPTION

%opts.XYStep = 6.055454452393343e-006
opts.XYStep = 3.0E-8;   % Transverse step (for findorbitx, findmxx,...)
opts.DPStep = 3.0E-6;   % Momentum step (for dispersion and chromaticity)
opts.OrbConvergence = 1.0E-12;  % Convergence of findorbitx
opts.OrbMaxIter = 20;           % Max number of iterations for findorbitx
end
