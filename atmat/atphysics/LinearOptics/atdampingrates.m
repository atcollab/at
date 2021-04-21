function [nu,chi]=atdampingrates(m66)
%ATDAMPINGRATES find tunes and damping rates from one map matrix with radiation
%
%[NU,CHI]=ATDAMPINGRATES(M66)
%
%note that in order to find the damping times, one needs the revolution
%time T0, then
%tau1 = T0/chi1, tau2 = T0/chi2, tau3 = T0/chi3

[~,vps]=amat(m66);
if length(vps) >= 3     % Invert rotation of longitudinale motion
    vps(3)=conj(vps(3));
end
nu=mod(angle(vps)/2/pi,1);
chi=-log(abs(vps));
end
