function [nu,chi]=atTunesAndDampingRatesFromMat(m66)
%find tunes and damping rates from one map matrix with radiation
%note that in order to find the damping times, one needs the revolution
%time T0, then
%tau1 = T0/chi1, tau2 = T0/chi2, tau3 = T0/chi3
%In the case that m66 is symplectic, the damping rates will be all 0.
%note that we have taken absolute values, assumeing all these quantities
%are positive.
%B. Nash 24/07/2014

aa=amat(m66);

Rmat=inv(aa)*m66*aa;

R1=Rmat([1 2],[1 2]);
R2=Rmat([3 4],[3 4]);
R3=Rmat([5 6],[5 6]);

%ev=eigs(m66);

ev1=eigs(R1);
evlog1=log(ev1(1));

ev2=eigs(R2);
evlog2=log(ev2(1));

ev3=eigs(R3);
evlog3=log(ev3(1));

%evlog123=log(ev([1,3,5]));

nu(1)=abs(imag(evlog1)/(2*pi));
chi(1)=abs(real(evlog1));

nu(2)=abs(imag(evlog2)/(2*pi));
chi(2)=abs(real(evlog2));

nu(3)=abs(imag(evlog3)/(2*pi));
chi(3)=abs(real(evlog3));