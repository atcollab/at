function sigz = sigz_zero(Vrf,circ,h,U0,E0,alpha,sigdelta)
% zero current bunch length

% Vrf is the RF voltage [V] (it may be a vector for multiple values)
% U0 is the energy loss around the ring [eV]
% h is the harmonic number
% alpha is the momentum compaction factor
% sigmadelta is the energy spread



phi=pi - asin(U0./Vrf);
nus= sqrt(-(Vrf/E0).*(h * alpha)/(2*pi) .* cos(phi));

sigz = circ*alpha/(2*pi.*nus)*sigdelta;
