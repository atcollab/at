function Q = Qval(Ib,Zn,Vrf,U0,E0,h,alpha,sigdelta)
% Qval gives the unitless Q parameter needed to compute the
% bunch lengthening factor due to the potential well effect

% Ib is the bunch current in mA
% Zn is the longitudinal broadband impedance
% Vrf is the RF voltage in MV
% h is the harmonic number
% alpha is the momentum compaction factor
% sigmadelta is the energy spread
% nus is the synchrotron tune

phi=pi - asin(U0/Vrf);
nus= sqrt(-(Vrf/E0)*(h * alpha)/(2*pi) * cos(phi));

Delta = -(2*pi*Ib*Zn)/(Vrf*h*cos(phi)*(alpha*sigdelta/nus)^3);
Q=Delta/(4*sqrt(pi));