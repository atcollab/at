function Nus = nus (VrfMV, alpha, U0MeV, E0MeV, h)

% this function return the synchrotron tune
% input:
% VrfMV is the RF voltage in MV
% alpha is the momentum compaction factor
% U0MeV is the energy lost per turn in MeV
% E0MeV is the beam energy in MeV
% h is the harmonic number

Phi=pi - asin(U0MeV./VrfMV);

Nus= sqrt(-(VrfMV/E0MeV)*(h * alpha)/(2*pi) .* cos(Phi));
