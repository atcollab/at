function rp = ringpara(THERING,varargin)
%RINGPARA Calculates various ring parameters
%(1) The calculation of emittance, mcf, momentum spread, bunch length, damping time, etc 
%is more accurate than atsummary.m because detailed
%calculation of dispersion function and curly H function inside dipoles is performed. 
%(2) calculate contribution of dispersion to vertical emittance.
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

if nargin==0
    global THERING; %#ok<TLEV>
end
Cq = 3.8319E-13; 
a = findcells(THERING,'Energy');
if isempty(a)
   gamma = 3000/PhysConstant.electron_mass_energy_equivalent_in_MeV.value;
else
   gamma = THERING{a(1)}.Energy/(PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6); 
end

dpindex = findcells(THERING,'BendingAngle');
[tw,tune,chrom] = twissring(THERING,0,dpindex,'chrom',0.00001);
beta = cat(1, tw.beta);
alpha = cat(1, tw.alpha);
mu = cat(1, tw.mu);
disper = cat(1, tw.Dispersion);
Dx  = disper(1:4:end);
Dxp = disper(2:4:end);
Dy  = disper(3:4:end);
Dyp = disper(4:4:end);

[tmptw,tune,chrom] = twissring(THERING,0,1:length(THERING),'chrom',0.00001);

lendp = getcellstruct(THERING,'Length',dpindex); %bending magnet length
lendp(lendp==0)=1;
theta = getcellstruct(THERING,'BendingAngle',dpindex); %bending angle
rho = lendp./theta;%THERING{dpindex(1)}.Length/THERING{dpindex(1)}.BendingAngle;

I1 = 0;
I2 = 0;
I3 = 0;
I4 = 0;
I5 = 0;

len = length(dpindex);
curHavg1 = 1:len;
for ii=1:len
  if theta(ii) ~= 0.0
      K = 0;
      Kk = 0;
      Kp = 0;
      if isfield(THERING{dpindex(ii)},'K')
          Kk = THERING{dpindex(ii)}.K;
      end
      if isfield(THERING{dpindex(ii)},'PolynomB')
          Kp = THERING{dpindex(ii)}.PolynomB(2);
      end
      if Kk~=Kp && (Kk~=0 && Kp~=0)
          Kk=0;
          warning('Values in K and PolynomB(2) are different and both not zero. Using PolynomB(2).'); 
      end
      Ks=[Kk,Kp];
      [~,i]=max(abs(Ks));
      K=Ks(i);
      
    th1 = THERING{dpindex(ii)}.EntranceAngle;
    th2 = THERING{dpindex(ii)}.ExitAngle;
    [dI1,dI2,dI3,dI4,dI5,curHavg1(ii)] = calcRadInt(rho(ii),theta(ii), ...
         alpha(ii,1),beta(ii,1),Dx(ii),Dxp(ii),K,th1,th2);
    I1 = I1 + dI1;
    I2 = I2 + dI2;
    I3 = I3 + dI3;
    I4 = I4 + dI4;
    I5 = I5 + dI5;
  end
end
% curHavg = sum(curHavg1.*lendp./abs(rho))/sum(lendp);
% %emittx =  Cq*gamma^2*curHavg/Jx/rho*1e9; %nm-rad
% emittx =  Cq*gamma^2*curHavg/Jx*1e9; %nm-rad
R = findspos(THERING, length(THERING)+1)/2/pi;
alphac = I1/2/pi/R;
U0 = 14.085*(gamma*PhysConstant.electron_mass_energy_equivalent_in_MeV.value/1000)^4*I2*1000.; %eV
if nargin>=2
    fprintf('dipole radiation loss:  %4.5f keV\n', U0/1000.);
    U0 = varargin{1}*1e6; %convert MeV to eV 
end
sigma_E = gamma*sqrt(Cq*I3/(2*I2+I4));
Jx = 1-I4/I2;
Jy = 1.00;
Je = 2+I4/I2;
emittx = Cq*gamma^2*I5/(I2-I4);

% minimum emittance due to radiation 1/gamma cone (Handbook, Chao, eq23, pag 211)
spos=findspos(THERING,1:length(THERING)+1);
[tmptw,tune,chrom] = twissring(THERING,0,1:length(THERING),'chrom');
betaall = cat(1, tmptw.beta);
meanbetayovers=sum(betaall(:,2)'.*diff(spos))/spos(end);%mean(betaall(:,2));%
meaninvr2=mean((1./rho.^2));%sum((1./rhos.^2).*diff(spos)/spos(end));%
meaninvr3=mean((1./abs(rho.^3)));%sum((1./abs(rhos.^3)).*diff(spos)/spos(end));%
emitty_lim = Cq*meanbetayovers/2/Jy*meaninvr3/meaninvr2;


cspeed = PhysConstant.speed_of_light_in_vacuum.value; %m/s
T0 = 2*pi*R/cspeed;
alpha0 = U0/1.0e6/2/T0/(gamma*PhysConstant.electron_mass_energy_equivalent_in_MeV.value);
alphax = Jx*alpha0;  %horizontal damping rate, 1/s
alphay = Jy*alpha0;
alphaE = Je*alpha0;

rp.E0 = gamma*PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1E6;
rp.R = R;
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
%[nu,chi]=atTunesAndDampingRatesFromMat(findm66(atradon(THERING)));
try
    [nu,chi]=atTunesAndDampingRatesFromMat(findm66((THERING)));
catch exc
    warning('failed coupled damping times computation');
    nu=[NaN,NaN,NaN];
    chi=[NaN,NaN,NaN];
end

rp.coupleddampingtime=T0./chi;

rp.dampingJ = [Jx,Jy,Je];

%rp.tw = tw;
%rp.tmptw = tmptw;
rp.tunes = tune;
rp.chroms = chrom;
rp.etac = 1/gamma^2-alphac;

cavind = findcells(THERING,'HarmNumber');
if ~isempty(cavind)
    freq_rf = THERING{cavind(1)}.Frequency;
    harm = THERING{cavind(1)}.HarmNumber;
    Vrf = 0;
    for ii=1:length(cavind)
        Vrf = Vrf+THERING{cavind(ii)}.Voltage;
    end
else
    % Default
    fprintf('warning: SPEAR3 rf parameters are assume\n');
    freq_rf = 476.314e6;
    harm = 372;
    Vrf = 3.2e6;
end

phi_s = pi-asin(rp.U0/Vrf);
nus = sqrt(harm*Vrf*abs(rp.etac*cos(phi_s))/2/pi/rp.E0);
rp.nus = nus;
rp.phi_s = phi_s;
rp.harm = harm;
rp.bunchlength = rp.sigma_E*rp.harm*abs(rp.etac)/rp.nus/2/pi/freq_rf*cspeed; % rms bunchlength in meter
delta_max = sqrt(2*U0/pi/alphac/harm/rp.E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));
rp.delta_max = delta_max;

%calculate vertical emittance
%1. contribution of vertical dispersion
curVavg1 = 1./beta(:,2).*(Dy.^2 + (beta(:,2).*Dyp + alpha(:,2).*Dy).^2);
curVavg = sum(curVavg1.*lendp./abs(rho))/sum(lendp);
emitty_d = Cq*gamma^2*curVavg/Jy; %m-rad

% %2. contribution of linear coupling resonance
% [G,Delta] = calc_lcG(THERING); 
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


function [dI1,dI2,dI3,dI4,dI5,curHavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1,th1,th2)
%[dI1,dI2,dI3,dI4,dI5,curHavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1)
%calcuate the contribution to the radiation integrals of a dipole. 
%rho, theta, radius and angle of the dipole
%a0, b0, are horizontal alpha and beta at the entrance of the dipole, 
%D0, D0p are dispersion at the entrace of the dipole
%K1, the gradient parameter in AT's convention, i.e., positive for
%horizontal focusing, K1=0 by default
%th1, th2, the entrance and exit angle, respectively, th1=th2=theta/2 by
%default. 
%

if nargin==6
   K1=0; 
end
if nargin<9
   th1 = 0; %theta/2.0;
   th2 = 0; %theta/2.0;
end

M21 = tan(th1)/rho;
D0p = M21*D0+D0p;
a0 = -M21*b0+a0;

N = 100;
th = (0:N)/N*theta;
len = length(th);
Dx = zeros(len,1); Dxp = zeros(len,1); curHavg1 = zeros(len,1);
for ii=1:len
       [Dx(ii), Dxp(ii)] = calcdisp(rho, th(ii), D0, D0p, K1);
       [ax, bx] = calctwiss(rho, th(ii), a0, b0, K1);
       curHavg1(ii) = (Dx(ii)^2+(ax*Dx(ii)+bx*Dxp(ii))^2)/bx;
end
curHavg = ((curHavg1(1)+curHavg1(end))/2.0 + sum(curHavg1(2:end-1)))/(length(th)-1);

dI1 = ((Dx(1) + Dx(end))/2.0 + sum(Dx(2:end-1)))*theta/N;
dI2 = abs(theta/rho);
dI3 = abs(theta/rho^2);
dI4 = (1/rho^2 + 2*K1)*dI1  - (Dx(1)/rho^2*tan(th1) + Dx(end)/rho^2*tan(th2));
dI5 = curHavg*abs(theta/rho^2);

function [Dx, Dxp] = calcdisp(rho, theta, D0, D0p, K1)
%calcualte dispersion function inside a combined-function dipole
s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Dx =  D0*cos(sqK*s) + D0p/sqK*sin(sqK*s)+(1-cos(sqK*s))/rho/sqK^2;
    Dxp = -D0*sqK*sin(sqK*s)+D0p*cos(sqK*s)+sin(sqK*s)/rho/sqK;
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Dx =  D0*cosh(sqK*s) + D0p/sqK*sinh(sqK*s)+(-1+cosh(sqK*s))/rho/sqK^2;
    Dxp = D0*sqK*sinh(sqK*s)+D0p*cosh(sqK*s)+sinh(sqK*s)/rho/sqK;

end

function [ax, bx] = calctwiss(rho, theta, a0, b0, K1)
%calculate twiss function inside a combined-function dipole manget
Mx = calcMx(rho, K1,theta);
g0 = (1+a0^2)/b0;
twx2 = [Mx(1,1)^2, -2*Mx(1,1)*Mx(1,2), Mx(1,2)^2; 
    -Mx(1,1)*Mx(2,1), Mx(1,1)*Mx(2,2)+Mx(1,2)*Mx(2,1),-Mx(1,2)*Mx(2,2);
    Mx(2,1)^2, -2*Mx(2,1)*Mx(2,2),Mx(2,2)^2] * [b0, a0, g0]';
ax = twx2(2);
bx = twx2(1);

function Mx = calcMx(rho,K1,theta)
s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Mx = [cos(sqK*s), sin(sqK*s)/sqK; -sqK*sin(sqK*s), cos(sqK*s)];
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Mx = [cosh(sqK*s), sinh(sqK*s)/sqK; sqK*sinh(sqK*s), cosh(sqK*s)];
end

