function rp = ringpara(varargin)
%RINGPARA Calculates various ring parameters
%
% rp = ringpara, use global THERING
% rp = ringpara(THERING)
% rp = ringpara(THERING,U0), supply total radiation loss in MeV
%
%  INPUTS
%  1. THERING - AT structure
%  2. DP - Energy offset
%
%  OUPUTS
%  1. RP - Structure with ring parameters
%
%  See also atx atsummary

%
%%Written by Xiaobiao Huang
%created on 12/17/2007
%Part of this code was modified from atsummary.m
%
%Modified by Peace Chang (check if theta(ii) ~= 0.0)
%Modified by S.Liuzzo and B.Nash (Dipole gradient may be in PolynomB(2),
%also coupled damping added) 7/24/2014
%
%Modified by A.Mash'al, Iranian Light Source Facility, 2018-07-30
%radiation effects of IDs added 

%Modified by Laurent Farvacque, 2020-08-27
%radiation effects removed until repaired...
%Analytical computation of radiation integrals

global THERING

e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
hbar=PhysConstant.Planck_constant_over2pi_times_c_in_MeV_fm.value;
cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s

Cq=55/32/sqrt(3)*hbar/e_mass*1E-9;   % m

[ring,Ux]=getargs(varargin,{THERING},[]);
try
    % Set values for 1 superperiod
    [energy,periods,Vrf0,harm0,U00]=atenergy(ring);
    Vrf=Vrf0/periods;
    harm=harm0/periods;
    U0=U00/periods;
catch ME
    warning(ME.identifier,'%s',ME.message);
    warning(ME.identifier,'Energy set to 3 GeV');
    energy=3E9;
    harm = 372;
    Vrf = 3.2e6;
end
if ~isempty(Ux)
    U0 = Ux*1e6; %convert MeV to eV 
    fprintf('dipole radiation loss:  %4.5f keV\n', U0/1000.);
end

gamma = energy/e_mass;
[lindata,tune,chrom]=atlinopt(ring,0.0,1:length(ring)+1);
[I1d,I2d,I3d,I4d,I5d,I6,Iv] = DipoleRadiation(ring,lindata);
[I1w,I2w,I3w,I4w,I5w] = WigglerRadiation(ring,lindata);
I1=I1d+I1w;
I2=I2d+I2w;
I3=I3d+I3w;
I4=I4d+I4w;
I5=I5d+I5w;

spos=findspos(ring,1:length(ring)+1);
circ=spos(end);

alphac = I1/circ;
sigma_E = gamma*sqrt(Cq*I3/(2*I2+I4));
Jx = 1-I4/I2;
Jy = 1.00;
Je = 2+I4/I2;
emittx = Cq*gamma^2*I5/(I2-I4);

% minimum emittance due to radiation 1/gamma cone (Handbook, Chao, eq23, pag 211)
emitty_lim = 13/55*Cq/Jy/I2*Iv;

T0 = circ/cspeed;
alpha0 = U0/2/T0/energy;
alphax = Jx*alpha0;  %horizontal damping rate, 1/s
alphay = Jy*alpha0;
alphaE = Je*alpha0;

rp.E0 = energy;
rp.R = circ/2/pi;
rp.alphac = alphac;
rp.U0 = U0; %eV
rp.sigma_E = sigma_E;
rp.emittx = emittx;
rp.emitty_lim= emitty_lim;
rp.T0 = T0;
rp.integrals = [I1,I2,I3,I4,I5];
rp.dampingalpha = [alphax, alphay, alphaE];
rp.dampingtime = 1./[alphax, alphay, alphaE];


%compute coupled damping times
%[nu,chi]=atTunesAndDampingRatesFromMat(findm66(atradon(ring)));
try
    [~,chi]=atTunesAndDampingRatesFromMat(findm66((ring)));
catch
    warning('failed coupled damping times computation');
    chi=[NaN,NaN,NaN];
end

rp.coupleddampingtime=T0./chi;

rp.dampingJ = [Jx,Jy,Je];

rp.tunes = tune;
rp.chroms = chrom;
rp.etac = 1/gamma^2-alphac;

cavind = findcells(ring,'HarmNumber');
if ~isempty(cavind)
    freq_rf = ring{cavind(1)}.Frequency;
else
    % Default
    fprintf('warning: SPEAR3 rf parameters are assume\n');
    freq_rf = 476.314e6;
end

phi_s = pi-asin(rp.U0/Vrf);
nus = sqrt(harm*Vrf*abs(rp.etac*cos(phi_s))/2/pi/rp.E0);
rp.nus = nus;
rp.phi_s = phi_s;
rp.harm = harm;
rp.bunchlength = rp.sigma_E*harm*abs(rp.etac)/nus/2/pi/freq_rf*cspeed; % rms bunchlength in meter
delta_max = sqrt(2*U0/pi/alphac/harm/rp.E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));
rp.delta_max = delta_max;

%calculate vertical emittance
%1. contribution of vertical dispersion
dipindex=atgetcells(ring,'BendingAngle');
lendp=atgetfieldvalues(ring(dipindex),'Length');
thedp=atgetfieldvalues(ring(dipindex),'BendingAngle');
beta=cat(1,lindata(dipindex).beta);
alpha=cat(1,lindata(dipindex).alpha);
disp=cat(2,lindata(dipindex).Dispersion)';
Dy=disp(:,3);
Dyp=disp(:,4);
curVavg1 = 1./beta(:,2).*(Dy.^2 + (beta(:,2).*Dyp + alpha(:,2).*Dy).^2);
curVavg = sum(curVavg1.*abs(thedp))/sum(lendp);
emitty_d = Cq*gamma^2*curVavg/Jy; %m-rad

% %2. contribution of linear coupling resonance
% [G,Delta] = calc_lcG(ring); 
% %emitty_c = emittx*abs(G)^2/(Delta^2+abs(G)^2);
% emitty_c = emittx*abs(G)^2/Delta^2/2.0; 
% rp.emitty_c = emitty_c;

rp.emitty_d = emitty_d;
% rp.emitty = emitty_d + emitty_c;

if nargout == 0
    fprintf('\n');
    fprintf('   ******** AT Ring Parmater Summary ********\n');
    fprintf('   Energy: \t\t\t%4.5f [GeV]\n', rp.E0/1E9);
    fprintf('   Circumference: \t\t%4.5f [m]\n', rp.R*2*pi);
    fprintf('   Revolution time: \t\t%4.5f [ns] (%4.5f [MHz]) \n', rp.T0*1e9,1./rp.T0*1e-6);
    fprintf('   Betatron tune H: \t\t%4.5f (%4.5f [kHz])\n', rp.tunes(1),(rp.tunes(1)-floor(rp.tunes(1)))/rp.T0*1e-3);
    fprintf('                 V: \t\t%4.5f (%4.5f [kHz])\n', rp.tunes(2),(rp.tunes(2)-floor(rp.tunes(2)))/rp.T0*1e-3);
    fprintf('   Momentum Compaction Factor: \t%4.5f\n', rp.alphac);
    fprintf('   Chromaticity H: \t\t%+4.5f\n', rp.chroms(1));
    fprintf('                V: \t\t%+4.5f\n', rp.chroms(2));
    fprintf('   Synchrotron Integral 1: \t%4.5f [m]\n', rp.integrals(1));
    fprintf('                        2: \t%4.5f [m^-1]\n', rp.integrals(2));
    fprintf('                        3: \t%4.5f [m^-2]\n', rp.integrals(3));
    fprintf('                        4: \t%4.5f [m^-1]\n', rp.integrals(4));
    fprintf('                        5: \t%4.5f [m^-1]\n', rp.integrals(5));
    fprintf('   Damping Partition H: \t%4.5f\n', rp.dampingJ(1));
    fprintf('                     V: \t%4.5f\n', rp.dampingJ(2));
    fprintf('                     E: \t%4.5f\n', rp.dampingJ(3));
    fprintf('   Radiation Loss: \t\t%4.5f [keV]\n', rp.U0/1000.);
    fprintf('   Natural Energy Spread: \t%4.5e\n', rp.sigma_E);
    fprintf('   Natural Emittance: \t\t%4.5e [nm]\n', rp.emittx*1e9);
    fprintf('   Radiation Damping H: \t%4.5f [ms]\n', rp.dampingtime(1)*1e3);
    fprintf('                     V: \t%4.5f [ms]\n', rp.dampingtime(2)*1e3);
    fprintf('                     E: \t%4.5f [ms]\n', rp.dampingtime(3)*1e3);
    fprintf('   Slip factor : \t%4.5f\n', rp.etac);
    fprintf('\n');
    fprintf('   Assuming cavities Voltage: %4.5f [kV]\n', Vrf/1e3);
    fprintf('                   Frequency: %4.5f [MHz]\n', freq_rf/1e6);
    fprintf('             Harmonic Number: %5d\n', rp.harm);
    fprintf('   Synchronous Phase:  %4.5f [rad] (%4.5f [deg])\n', rp.phi_s, rp.phi_s*180/pi);
    fprintf('   Linear Energy Acceptance:  %4.5f %%\n', rp.delta_max*100);
    fprintf('   Synchrotron Tune:   %4.5f (%4.5f kHz or %4.2f turns) \n', rp.nus, rp.nus/rp.T0*1e-3, 1/rp.nus);
    fprintf('   Bunch Length:       %4.5f [mm], %4.5f [ps]\n', rp.bunchlength*1e3, rp.bunchlength/cspeed*1e12);
    fprintf('\n');
%     fprintf('   Vertical Emittance:  %4.5f [nm]\n', rp.emitty*1e9);
%     fprintf('   Emitty from Dy:  %4.5f [nm], from linear coupling: %4.5f\n', rp.emitty_d*1e9,rp.emitty_c*1e9);
    fprintf('   Emitty from Dy:  %4.5f [nm]\n', rp.emitty_d*1e9);
    fprintf('   Emitty 1/gamma cone limit:  %4.5f [pm]\n', rp.emitty_lim*1e12);
end
