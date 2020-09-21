function varargout = atsummary(varargin)
%ATSUMMARY - Prints out the paramters of the current AT lattice
%  The parameters that come after the Synchrotron Integrals are
%  parameters that depend on the Integrals themselves. The equations to
%  calculate them were taken from [1].
%
%  [1] Alexander Wu Chao and Maury Tigner, Handbook of Accelerator Physics
%  and Engineering (World Scientific, Singapore, 1998), pp. 183-187. (or
%  187-190 in ed. 2)

%
%  Written by Eugene Tan
%  Revised by Laurent S. Nadolski

global THERING GLOBVAL

e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
e_radius=PhysConstant.classical_electron_radius.value;                  % m
hbar=PhysConstant.Planck_constant_over2pi_times_c_in_MeV_fm.value;
cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s

Cq=55/32/sqrt(3)*hbar/e_mass*1E-9;                                      % m
Cgamma=4.0E27*pi*e_radius/3/e_mass^3;                                   % [m/GeV^3]

[nodisp,varargs]=getflag(varargin,'NoDisplay');
[disp,varargs]=getflag(varargs,'Display');
ring=getargs(varargs,{THERING});
DisplayFlag=(~nodisp) || disp;

energy=atenergy(ring);
params   = atgetcells(ring,'Class','RingParam');
if any(params)
    LatticeName=ring{find(params,1)}.FamName;
elseif isfield(GLOBVAL,'LatticeFile')
    LatticeName=GLOBVAL.LatticeFile;
else
    LatticeName='';
end

cavities = atgetcells(ring,'Frequency');
if any(cavities)
    freq = ring{find(cavities,1)}.Frequency;
    v_cav=sum(atgetfieldvalues(ring(cavities),'Voltage'));
else
    % Default
    freq = 352.202e6;
    v_cav = 3e6;
end

% Structure to store info
smm.e0            = energy*1e-9;
smm.circumference = findspos(ring, length(ring)+1);
smm.revTime       = smm.circumference / cspeed;
smm.revFreq       = cspeed / smm.circumference;
smm.gamma         = 1.e9 * smm.e0 / e_mass;
smm.beta          = sqrt(1 - 1/smm.gamma);

[TD, smm.tunes, smm.chromaticity] = twissring(ring, 0, 1:length(ring)+1, 'chrom', 1e-8);
smm.compactionFactor = mcf(ring);

% For calculating the synchrotron integrals
[I1,I2,I3,I4,I5,I6,~] = DipoleRadiation(ring,TD);
smm.integrals=[I1,I2,I3,I4,I5,I6];

% Damping numbers
% Use Robinson's Theorem
smm.damping(1) = 1 - smm.integrals(4)/smm.integrals(2);
smm.damping(2) = 1;
smm.damping(3) = 2 + smm.integrals(4)/smm.integrals(2);

smm.radiation           = Cgamma/2/pi*smm.e0.^4*smm.integrals(2);    % GeV
smm.naturalEnergySpread = smm.gamma*sqrt(Cq*smm.integrals(3)/(2*smm.integrals(2) + smm.integrals(4)));
smm.naturalEmittance    = Cq*smm.gamma.^2*smm.integrals(5)/(smm.integrals(2)-smm.integrals(4));

% Damping times
cf=smm.radiation/2/smm.revTime/smm.e0;
smm.radiationDamping(1) = 1/(cf*smm.damping(1));
smm.radiationDamping(2) = 1/(cf*smm.damping(2));
smm.radiationDamping(3) = 1/(cf*smm.damping(3));

% Slip factor
smm.etac = smm.gamma^(-2) - smm.compactionFactor;

smm.harmon = smm.circumference/(PhysConstant.speed_of_light_in_vacuum.value/freq); % Assuming 499.654MHz RF %
smm.overvoltage = v_cav/(smm.radiation*1e9); 
% Assuming the harmon and overvoltage above.
% references:  H. Winick, "Synchrotron Radiation Sources: A Primer",
% World Scientific Publishing, Singapore, pp92-95. (1995)
% Wiedemann, pp290,350. Chao, pp189.
smm.syncphase = pi - asin(1/smm.overvoltage);
smm.energyacceptance = sqrt(v_cav*sin(smm.syncphase)*2*(sqrt(smm.overvoltage^2-1) - acos(1/smm.overvoltage))/(pi*smm.harmon*abs(smm.etac)*smm.e0*1e9));
smm.synctune = sqrt((smm.etac*smm.harmon*v_cav*cos(smm.syncphase))/(2*pi*smm.e0*1e9));
smm.bunchlength = smm.beta*cspeed*abs(smm.etac)*smm.naturalEnergySpread/(smm.synctune*smm.revFreq*2*pi);

% optics
% [bx by] = modelbeta;
% [ax ay] = modelalpha;
% [etax etay] = modeleta;
% [etaprimex etaprimey] = modeletaprime;

bx = TD(1).beta(1);
by = TD(1).beta(2);
ax = TD(1).alpha(1);
ay = TD(1).alpha(2);
etax = TD(1).Dispersion(1);
etay = TD(1).Dispersion(3);
etaprimex = TD(1).Dispersion(2);
etaprimey = TD(1).Dispersion(4);

if DisplayFlag
    SeparatorString = '   ******************************************************************\n';
    fprintf('\n');
    fprintf('   *************  Summary for ''%s'' ************\n', LatticeName);
    fprintf('   Energy: \t\t\t% 4.5f [GeV]\n', smm.e0);
    fprintf('   Gamma: \t\t\t% 4.5f \n', smm.gamma);
    fprintf('   Circumference: \t\t% 4.5f [m]\n', smm.circumference);
    fprintf('   Revolution time: \t\t% 4.5f [ns] (%4.5f [MHz]) \n', smm.revTime*1e9,smm.revFreq*1e-6);
    fprintf('   Betatron tune H: \t\t% 4.5f (%4.5f [kHz])\n', smm.tunes(1),smm.tunes(1)/smm.revTime*1e-3);
    fprintf('                 V: \t\t% 4.5f (%4.5f [kHz])\n', smm.tunes(2),smm.tunes(2)/smm.revTime*1e-3);
    fprintf('   Momentum Compaction Factor: \t% 4.5e\n', smm.compactionFactor);
    fprintf('   Chromaticity H: \t\t%+4.5f\n', smm.chromaticity(1));
    fprintf('                V: \t\t%+4.5f\n', smm.chromaticity(2));
    fprintf(SeparatorString);
    fprintf('   Synchrotron Integral 1: \t% 4.5e [m]\n', smm.integrals(1));
    fprintf('                        2: \t% 4.5e [m^-1]\n', smm.integrals(2));
    fprintf('                        3: \t% 4.5e [m^-2]\n', smm.integrals(3));
    fprintf('                        4: \t% 4.5e [m^-1]\n', smm.integrals(4));
    fprintf('                        5: \t% 4.5e [m^-1]\n', smm.integrals(5));
    fprintf('                        6: \t% 4.5e [m^-1]\n', smm.integrals(6));
    fprintf('   Damping Partition H: \t% 4.5f\n', smm.damping(1));
    fprintf('                     V: \t% 4.5f\n', smm.damping(2));
    fprintf('                     E: \t% 4.5f\n', smm.damping(3));
    fprintf('   Radiation Loss: \t\t% 4.5f [keV]\n', smm.radiation*1e6);
    fprintf('   Natural Energy Spread: \t% 4.5e\n', smm.naturalEnergySpread);
    fprintf('   Natural Emittance: \t\t% 4.5e [mrad]\n', smm.naturalEmittance);
    fprintf('   Radiation Damping H: \t% 4.5f [ms] or %4.2f turns\n', smm.radiationDamping(1)*1e3, smm.radiationDamping(1)/smm.revTime);
    fprintf('                     V: \t% 4.5f [ms] or %4.2f turns\n', smm.radiationDamping(2)*1e3, smm.radiationDamping(2)/smm.revTime);
    fprintf('                     E: \t% 4.5f [ms] or %4.2f turns\n', smm.radiationDamping(3)*1e3, smm.radiationDamping(3)/smm.revTime);
    fprintf('   Slip factor: \t\t%4.5e\n', smm.etac);
    fprintf('   Momentum compaction factor: \t %4.5e (%4.5e)\n',  smm.integrals(1)/smm.circumference, smm.compactionFactor);
    fprintf(SeparatorString);
    fprintf('   Assuming cavities Voltage: \t% 4.5f [kV]\n', v_cav/1e3);
    fprintf('                   Frequency: \t% 4.5f [MHz]\n', freq/1e6);
    fprintf('             Harmonic Number: \t% 4.0f\n', smm.harmon);
    fprintf('   Overvoltage factor: \t\t% 4.5f\n', smm.overvoltage);
    fprintf('   Synchronous Phase:  \t\t% 4.5f [rad] (%4.5f [deg])\n', smm.syncphase, smm.syncphase*180/pi);
    fprintf('   Linear Energy Acceptance:  \t% 4.3f %%\n', smm.energyacceptance*100);
    fprintf('   Synchrotron Tune:   \t\t% 4.5f (%4.5f kHz or %4.2f turns) \n', smm.synctune, smm.synctune/smm.revTime*1e-3, 1/smm.synctune);
    fprintf('   Bunch Length:       \t\t% 4.5f [mm] (%4.5f ps)\n', smm.bunchlength*1e3, smm.bunchlength/3e8*1e12);
    fprintf(SeparatorString);
    fprintf('   Injection:\n');
    fprintf('   H: beta = %06.3f [m] alpha = %+04.1e eta = %+04.3f [m] eta'' = %+04.1e \n', bx(1), ax(1), etax(1), etaprimex(1));
    fprintf('   V: beta = %06.3f [m] alpha = %+04.1e eta = %+04.3f [m] eta'' = %+04.1e \n', by(1), ay(1), etay(1), etaprimey(1));
    fprintf('   ********** End of Summary for ''%s'' **********\n', LatticeName);
    fprintf('\n');
end

if nargout > 0
    varargout{1} = smm;
end
end
