function sum = atsummary
%ATSUMMARY - Prints out the paramters of the current AT lattice
%  The parameters that come after the Synchrotron Integrals are
%  parameters that depend on the Integrals themselves. The equations to
%  calculate them were taken from [1].
%
%  [1] Alexander Wu Chao and Maury Tigner, Handbook of Accelerator Physics
%  and Engineering (World Scientific, Singapore, 1998), pp. 183-187. (or
%  187-190 in ed. 2)
%
%  See also ringpara

%  Written by Eugene Tan
%  Revised by Laurent S. Nadolski


global THERING

% Structure to store info
sum.e0 = getenergy('Model');
sum.circumference = findspos(THERING, length(THERING)+1);
sum.revTime = sum.circumference / 2.99792458e8;
sum.revFreq = 2.99792458e8 / sum.circumference;
sum.gamma = sum.e0 / 0.51099906e-3;
sum.beta = sqrt(1 - 1/sum.gamma);
[TD, sum.tunes, sum.chromaticity] = twissring(THERING, 0, 1:length(THERING)+1, 'chrom', 1e-8);
sum.compactionFactor = mcf(THERING);

% For calculating the synchrotron integrals
temp  = cat(2,TD.Dispersion);
D_x   = temp(1,:)';
D_x_  = temp(2,:)';
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
gamma = (1 + alpha.^2) ./ beta;
circ  = TD(length(THERING)+1).SPos;

% Synchrotron integral calculation
sum.integrals = [0.0 0.0 0.0 0.0 0.0 0.0];

for i = 1:length(THERING),
    if isfield(THERING{i}, 'BendingAngle') && isfield(THERING{i}, 'EntranceAngle')
        rho = THERING{i}.Length/THERING{i}.BendingAngle;
        dispersion = 0.5*(D_x(i)+D_x(i+1));
        sum.integrals(1) = sum.integrals(1) + dispersion*THERING{i}.Length/rho;
        sum.integrals(2) = sum.integrals(2) + THERING{i}.Length/(rho^2);
        sum.integrals(3) = sum.integrals(3) + THERING{i}.Length/(rho^3);
        % For general wedge magnets
        sum.integrals(4) = sum.integrals(4) + ...
            D_x(i)*tan(THERING{i}.EntranceAngle)/rho^2 + ...
            (1 + 2*rho^2*THERING{i}.PolynomB(2))*(D_x(i)+D_x(i+1))*THERING{i}.Length/(2*rho^3) + ...
            D_x(i+1)*tan(THERING{i}.ExitAngle)/rho^2;
        %         sum.integrals(4) = sum.integrals(4) + 2*0.5*(D_x(i)+D_x(i+1))*THERING{i}.Length/rho^3;
        H1 = beta(i,1)*D_x_(i)*D_x_(i)+2*alpha(i)*D_x(i)*D_x_(i)+gamma(i)*D_x(i)*D_x(i);
        H0 = beta(i+1,1)*D_x_(i+1)*D_x_(i+1)+2*alpha(i+1)*D_x(i+1)*D_x_(i+1)+gamma(i+1)*D_x(i+1)*D_x(i+1);
        sum.integrals(5) = sum.integrals(5) + THERING{i}.Length*(H1+H0)*0.5/(rho^3);
        %         if H1+H0 < 0
        %             fprintf('%f %i %s\n', H1+H0, i, THERING{i}.FamName)
        %         end
        sum.integrals(6) = sum.integrals(6) + THERING{i}.PolynomB(2)^2*dispersion^2*THERING{i}.Length;
    end
end

% Damping numbers
% Use Robinson's Theorem
sum.damping(1) = 1 - sum.integrals(4)/sum.integrals(2);
sum.damping(2) = 1;
sum.damping(3) = 2 + sum.integrals(4)/sum.integrals(2);

sum.radiation = 8.846e-5*sum.e0.^4*sum.integrals(2)/(2*pi);
sum.naturalEnergySpread = sqrt(3.8319e-13*sum.gamma.^2*sum.integrals(3)/(2*sum.integrals(2) + sum.integrals(4)));
sum.naturalEmittance = 3.8319e-13*(sum.e0*1e3/0.510999).^2*sum.integrals(5)/(sum.damping(1)*sum.integrals(2));

% Damping times
sum.radiationDamping(1) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(1)/circ);
sum.radiationDamping(2) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(2)/circ);
sum.radiationDamping(3) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(3)/circ);

% Slip factor
sum.etac = sum.gamma^(-2) - sum.compactionFactor;

cavind = findcells(THERING,'HarmNumber');
if ~isempty(cavind)
    freq = THERING{cavind}.Frequency;
    v_cav = THERING{cavind}.Voltage;
else
    % Default
    freq = 352.202e6;
    v_cav = 3e6;
end
sum.harmon = sum.circumference/(2.99792458e8/freq); % Assuming 499.654MHz RF
sum.overvoltage = v_cav/(sum.radiation*1e9); % Assuming 3e6 volt cavities.
% Assuming the harmon and overvoltage above.
% references:  H. Winick, "Synchrotron Radiation Sources: A Primer",
% World Scientific Publishing, Singapore, pp92-95. (1995)
% Wiedemann, pp290,350. Chao, pp189.
sum.syncphase = pi - asin(1/sum.overvoltage);
sum.energyacceptance = sqrt(v_cav*sin(sum.syncphase)*2*(sqrt(sum.overvoltage^2-1) - acos(1/sum.overvoltage))/(pi*sum.harmon*abs(sum.etac)*sum.e0*1e9));
sum.synctune = sqrt((sum.etac*sum.harmon*v_cav*cos(sum.syncphase))/(2*pi*sum.e0*1e9));
sum.bunchlength = sum.beta*299792458*abs(sum.etac)*sum.naturalEnergySpread/(sum.synctune*sum.revFreq*2*pi);

if nargout == 0
    fprintf('\n');
    %fprintf('   ******** Summary for ''%s'' ********\n', GLOBVAL.LatticeFile);
    fprintf('   ******** AT Lattice Summary ********\n');
    fprintf('   Energy: \t\t\t%4.5f [GeV]\n', sum.e0);
    fprintf('   Gamma: \t\t\t%4.5f \n', sum.gamma);
    fprintf('   Circumference: \t\t%4.5f [m]\n', sum.circumference);
    fprintf('   Revolution time: \t\t%4.5f [ns] (%4.5f [MHz]) \n', sum.revTime*1e9,sum.revFreq*1e-6);
    fprintf('   Betatron tune H: \t\t%4.5f (%4.5f [kHz])\n', sum.tunes(1),sum.tunes(1)/sum.revTime*1e-3);
    fprintf('                 V: \t\t%4.5f (%4.5f [kHz])\n', sum.tunes(2),sum.tunes(2)/sum.revTime*1e-3);
    fprintf('   Momentum Compaction Factor: \t%4.5f\n', sum.compactionFactor);
    fprintf('   Chromaticity H: \t\t%+4.5f\n', sum.chromaticity(1));
    fprintf('                V: \t\t%+4.5f\n', sum.chromaticity(2));
    fprintf('   Synchrotron Integral 1: \t%4.5f [m]\n', sum.integrals(1));
    fprintf('                        2: \t%4.5f [m^-1]\n', sum.integrals(2));
    fprintf('                        3: \t%4.5f [m^-2]\n', sum.integrals(3));
    fprintf('                        4: \t%4.5f [m^-1]\n', sum.integrals(4));
    fprintf('                        5: \t%4.5f [m^-1]\n', sum.integrals(5));
    fprintf('                        6: \t%4.5f [m^-1]\n', sum.integrals(6));
    fprintf('   Damping Partition H: \t%4.5f\n', sum.damping(1));
    fprintf('                     V: \t%4.5f\n', sum.damping(2));
    fprintf('                     E: \t%4.5f\n', sum.damping(3));
    fprintf('   Radiation Loss: \t\t%4.5f [keV]\n', sum.radiation*1e6);
    fprintf('   Natural Energy Spread: \t%4.5e\n', sum.naturalEnergySpread);
    fprintf('   Natural Emittance: \t\t%4.5e [mrad]\n', sum.naturalEmittance);
    fprintf('   Radiation Damping H: \t%4.5f [ms]\n', sum.radiationDamping(1)*1e3);
    fprintf('                     V: \t%4.5f [ms]\n', sum.radiationDamping(2)*1e3);
    fprintf('                     E: \t%4.5f [ms]\n', sum.radiationDamping(3)*1e3);
    fprintf('   Slip factor : \t%4.5f\n', sum.etac);
    fprintf('\n');
    fprintf('   Assuming cavities Voltage: %4.5f [kV]\n', v_cav/1e3);
    fprintf('                   Frequency: %4.5f [MHz]\n', freq/1e6);
    fprintf('             Harmonic Number: %4.5f\n', sum.harmon);
    fprintf('   Overvoltage factor: %4.5f\n', sum.overvoltage);
    fprintf('   Synchronous Phase:  %4.5f [rad] (%4.5f [deg])\n', sum.syncphase, sum.syncphase*180/pi);
    fprintf('   Linear Energy Acceptance:  %4.5f %%\n', sum.energyacceptance*100);
    fprintf('   Synchrotron Tune:   %4.5f (%4.5f kHz or %4.2f turns) \n', sum.synctune, sum.synctune/sum.revTime*1e-3, 1/sum.synctune);
    fprintf('   Bunch Length:       %4.5f [mm]\n', sum.bunchlength*1e3);
end