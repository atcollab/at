function [I1,I2,I3,I4,I5,I6,Iv] = DipoleRadiation(ring,lindata)
%DIPOLERADIATION	Compute the radiation integrals in dipoles

% Analytical integration from:
%
% EVALUATION OF SYNCHROTRON RADIATION INTEGRALS
% R.H. Helm, M.J. Lee, P.L. Morton and M. Sands
% SLAC-PUB-1193, March 1973

[I1,I2,I3,I4,I5,I6,Iv] = ElementRadiation(ring,lindata,'UseQuadrupole', false);

