function delta_max_rf = RFacc(Vrf,U0,E0,h,alpha)
%RFACC Computes the RF acceptance with linear formula
%   delta_max_rf = RFacc(Vrf,U0,E0,h,alpha)
%
%   This function computes the RF acceptance
%   Vrf is the RF voltage in V
%   U0 is the energy loss per turn in eV
%   E0 is the energy of the beam in eV
%   h is the harmonic number
%   alpha is the momentum compaction factor
%
%  See also atRFacc

delta_max_rf = sqrt(2*U0/pi./alpha/h/E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));