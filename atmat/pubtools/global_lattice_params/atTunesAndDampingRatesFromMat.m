function [nu1,nu2,nu3,chi1,chi2,chi3]=atTunesAndDampingRatesFromMat(m66)
%find tunes and damping rates from one map matrix with radiation
%note that in order to find the damping times, one needs the revolution
%time T0, then
%tau1 = T0/chi1, tau2 = T0/chi2, tau3 = T0/chi3
%In the case that m66 is symplectic, the damping rates will be all 0.
%note that we have taken absolute values, assumeing all these quantities
%are positive.
%B. Nash 24/07/2014

ev=eigs(m66);

evlog123=log(ev([1,3,5]));

nu1=abs(imag(evlog123(1))/(2*pi));
chi1=abs(real(evlog123(1)));

nu2=abs(imag(evlog123(2))/(2*pi));
chi2=abs(real(evlog123(2)));

nu3=abs(imag(evlog123(3))/(2*pi));
chi3=abs(real(evlog123(3)));