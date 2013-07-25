function [tvac,sigmacoulomb,sigmabremsstrahlung]=...
    VacLifetimeResidualGass(ring,DEE,Zperc,P,b)
% Coloumb scattering and residual gas Bremsstrahlung Cross sections are
% computed to evaluate the lifetime contribution of the vacuum system
% specified by Zperc (vector of percentages, index=atomic Number),and P[nTor] the
% pressure profile along the ring
% b is the dimension of the vacuum chamber along the lattice (x,y).
% DEE is the Energy acceptance of the ring (longitudinal aperture)
%
% length(b)=length(DEE)=length(P)=length(ring)
%
% Assumed temperature of 293 K and interacion with diatomic gas
% 
% tvac returned in hours
% 
% created 20-5-2013

r0=2.8179403267e-15;%[m] classical electron radius
c=299792458;% m/s speed of ligth
me=0.510998910e6;% eV/c^2 electron mass
alpha=1/137.035999084; % fine structure constant
kboltz=1.3806488e-23 ; %m^2 kg s^-2 K^-1 Boltzman constant
T=293;%K

% particle energy
E0=ring{1}.Energy;% eV/c electron beam energy

l=atlinopt(ring,0,1:length(ring));
betx=arrayfun(@(x)x.beta(1),l)';
bety=arrayfun(@(x)x.beta(2),l)';

Z=1:length(Zperc);

sgcv=max(bety).*bety./b(:,2).^2;

sgch=max(betx).*betx./b(:,1).^2;

% test parameters, with this it should return 0.13 barn
% sgch=0;
% sgcv=30*15/0.03^2;
% E0=5e9;

sgc=sgch+sgcv; % this is CORRECT ????

ZavC=sum(Z.^2.*Zperc);% this is CORRECT ????

sigmacoulombLocal=2*pi.*r0.^2.*ZavC.*sgc.*(me/E0).^2;% sum contribution for various atoms

sigmacoulomb=mean(sigmacoulombLocal);

sigmacoulomb=sigmacoulomb*10^4*10^24;% cross section from m^2 to barn, 1 barn=10^-24 cm^2


% DEE is along the lattice. average of sigma along the lattice
Zav=sum(Z.*(Z+1).*log(183./Z.^(1/3)).*Zperc);

sigmabremsstrahlungLocal=16*alpha*r0^2/3.*Zav.*(log(1./DEE)-5/8);

sigmabremsstrahlung=mean(sigmabremsstrahlungLocal);

sigmabremsstrahlung=sigmabremsstrahlung*10^4*10^24;% cross section from m^2 to barn, 1 barn=10^-24 cm^2
size(sigmabremsstrahlungLocal)
size(sigmacoulombLocal)

sigma=sigmabremsstrahlungLocal*10^4*10^24+sigmacoulombLocal*10^4*10^24; %in barn

% n_mol=P./kboltz/T;%[m^-3] molar density of residual gas along the ring

invtvac=P.*sigma./(0.474*T);

tvac=length(sigma)/sum(invtvac);% h


