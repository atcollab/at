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

DisplayFlag = 1;

for i = length(varargin):-1:1
    if strcmpi(varargin(i),'NoDisplay')
        DisplayFlag = 0;
        varargin(i) = [];
    elseif strcmpi(varargin(i),'Display')
        DisplayFlag = 1;
        varargin(i) = [];
    end
end

if isempty(varargin)
    global THERING GLOBVAL;
else
    THERING = varargin{i};
    GLOBVAL.E0 = THERING{1}.Energy;
    GLOBVAL.LatticeFile = '';  
end

% Structure to store info
sum.e0            = GLOBVAL.E0*1e-9;
sum.circumference = findspos(THERING, length(THERING)+1);
sum.revTime       = sum.circumference / PhysConstant.speed_of_light_in_vacuum.value;
sum.revFreq       = PhysConstant.speed_of_light_in_vacuum.value / sum.circumference;
sum.gamma         = sum.e0 / PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e3; %0.51099906e-3; 
sum.beta          = sqrt(1 - 1/sum.gamma);

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
sum.integrals = zeros(1,6);

ii =0;
for i = 1:length(THERING)
    if isfield(THERING{i}, 'BendingAngle') && isfield(THERING{i}, 'EntranceAngle')
        ii = ii +1;
        rho = THERING{i}.Length/THERING{i}.BendingAngle;
        [dI1,dI2,dI3,dI4,dI5,curHavg1(ii), Dxavg(ii)] = ...
            calcRadInt(rho,THERING{i}.BendingAngle, ...
            alpha(i,1),beta(i,1),D_x(i),D_x_(i),...
            THERING{i}.K,THERING{i}.EntranceAngle,THERING{i}.ExitAngle);
        
        sum.integrals(1) = sum.integrals(1) + dI1;
        sum.integrals(2) = sum.integrals(2) + dI2;
        sum.integrals(3) = sum.integrals(3) + dI3;
        % For general wedge magnets
        sum.integrals(4) = sum.integrals(4) + dI4;
        sum.integrals(5) = sum.integrals(5) + dI5;
        %         sum.integrals(4) = sum.integrals(4) + 2*0.5*(D_x(i)+D_x(i+1))*THERING{i}.Length/rho^3;
        H1 = beta(i,1)*D_x_(i)*D_x_(i)+2*alpha(i)*D_x(i)*D_x_(i)+gamma(i)*D_x(i)*D_x(i);
        H0 = beta(i+1,1)*D_x_(i+1)*D_x_(i+1)+2*alpha(i+1)*D_x(i+1)*D_x_(i+1)+gamma(i+1)*D_x(i+1)*D_x(i+1);
        sum.integrals(6) = sum.integrals(6) + THERING{i}.PolynomB(2)^2*Dxavg(ii)^2*THERING{i}.Length;
    end
end

% Damping numbers
% Use Robinson's Theorem
sum.damping(1) = 1 - sum.integrals(4)/sum.integrals(2);
sum.damping(2) = 1;
sum.damping(3) = 2 + sum.integrals(4)/sum.integrals(2);

sum.radiation           = 8.846e-5*sum.e0.^4*sum.integrals(2)/(2*pi);
sum.naturalEnergySpread = sqrt(3.8319e-13*sum.gamma.^2*sum.integrals(3)/(2*sum.integrals(2) + sum.integrals(4)));
sum.naturalEmittance    = 3.8319e-13*(sum.e0*1e3/PhysConstant.electron_mass_energy_equivalent_in_MeV.value).^2*sum.integrals(5)/(sum.damping(1)*sum.integrals(2)); %% need to be replaced by constant ?

% Damping times
sum.radiationDamping(1) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(1)/circ);
sum.radiationDamping(2) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(2)/circ);
sum.radiationDamping(3) = 1/(2113.1*sum.e0.^3*sum.integrals(2)*sum.damping(3)/circ);

% Slip factor
sum.etac = sum.gamma^(-2) - sum.compactionFactor;

cavind = findcells(THERING,'HarmNumber');
if ~isempty(cavind)
    freq = THERING{cavind(:,1)}.Frequency;
    v_cav = 0;
    for i = 1:length(cavind)
        v_cav = v_cav + THERING{cavind(:,1)}.Voltage;
    end
else
    % Default
    freq = 352.202e6;
    v_cav = 3e6;
end
sum.harmon = sum.circumference/(PhysConstant.speed_of_light_in_vacuum.value/freq); % Assuming 499.654MHz RF %
sum.overvoltage = v_cav/(sum.radiation*1e9); 
% Assuming the harmon and overvoltage above.
% references:  H. Winick, "Synchrotron Radiation Sources: A Primer",
% World Scientific Publishing, Singapore, pp92-95. (1995)
% Wiedemann, pp290,350. Chao, pp189.
sum.syncphase = pi - asin(1/sum.overvoltage);
sum.energyacceptance = sqrt(v_cav*sin(sum.syncphase)*2*(sqrt(sum.overvoltage^2-1) - acos(1/sum.overvoltage))/(pi*sum.harmon*abs(sum.etac)*sum.e0*1e9));
sum.synctune = sqrt((sum.etac*sum.harmon*v_cav*cos(sum.syncphase))/(2*pi*sum.e0*1e9));
sum.bunchlength = sum.beta*PhysConstant.speed_of_light_in_vacuum.value*abs(sum.etac)*sum.naturalEnergySpread/(sum.synctune*sum.revFreq*2*pi);

% optics
% [bx by] = modelbeta;
% [ax ay] = modelalpha;
% [etax etay] = modeleta;
% [etaprimex etaprimey] = modeletaprime;

[LinData,~, ~] = atlinopt(THERING,0,1);
bx = LinData.beta(1);
by = LinData.beta(2);
ax = LinData.alpha(1);
ay = LinData.alpha(2);
etax = LinData.Dispersion(1);
etay = LinData.Dispersion(3);
etaprimex = LinData.Dispersion(2);
etaprimey = LinData.Dispersion(4);

if DisplayFlag
    SeparatorString = '   ******************************************************************\n';
    fprintf('\n');
    fprintf('   *************  Summary for ''%s'' ************\n', GLOBVAL.LatticeFile);
    fprintf('   Energy: \t\t\t% 4.5f [GeV]\n', sum.e0);
    fprintf('   Gamma: \t\t\t% 4.5f \n', sum.gamma);
    fprintf('   Circumference: \t\t% 4.5f [m]\n', sum.circumference);
    fprintf('   Revolution time: \t\t% 4.5f [ns] (%4.5f [MHz]) \n', sum.revTime*1e9,sum.revFreq*1e-6);
    fprintf('   Betatron tune H: \t\t% 4.5f (%4.5f [kHz])\n', sum.tunes(1),sum.tunes(1)/sum.revTime*1e-3);
    fprintf('                 V: \t\t% 4.5f (%4.5f [kHz])\n', sum.tunes(2),sum.tunes(2)/sum.revTime*1e-3);
    fprintf('   Momentum Compaction Factor: \t% 4.5e\n', sum.compactionFactor);
    fprintf('   Chromaticity H: \t\t%+4.5f\n', sum.chromaticity(1));
    fprintf('                V: \t\t%+4.5f\n', sum.chromaticity(2));
    fprintf(SeparatorString);
    fprintf('   Synchrotron Integral 1: \t% 4.5e [m]\n', sum.integrals(1));
    fprintf('                        2: \t% 4.5e [m^-1]\n', sum.integrals(2));
    fprintf('                        3: \t% 4.5e [m^-2]\n', sum.integrals(3));
    fprintf('                        4: \t% 4.5e [m^-1]\n', sum.integrals(4));
    fprintf('                        5: \t% 4.5e [m^-1]\n', sum.integrals(5));
    fprintf('                        6: \t% 4.5e [m^-1]\n', sum.integrals(6));
    fprintf('   Damping Partition H: \t% 4.5f\n', sum.damping(1));
    fprintf('                     V: \t% 4.5f\n', sum.damping(2));
    fprintf('                     E: \t% 4.5f\n', sum.damping(3));
    fprintf('   Radiation Loss: \t\t% 4.5f [keV]\n', sum.radiation*1e6);
    fprintf('   Natural Energy Spread: \t% 4.5e\n', sum.naturalEnergySpread);
    fprintf('   Natural Emittance: \t\t% 4.5e [mrad]\n', sum.naturalEmittance);
    fprintf('   Radiation Damping H: \t% 4.5f [ms] or %4.2f turns\n', sum.radiationDamping(1)*1e3, sum.radiationDamping(1)/sum.revTime);
    fprintf('                     V: \t% 4.5f [ms] or %4.2f turns\n', sum.radiationDamping(2)*1e3, sum.radiationDamping(2)/sum.revTime);
    fprintf('                     E: \t% 4.5f [ms] or %4.2f turns\n', sum.radiationDamping(3)*1e3, sum.radiationDamping(3)/sum.revTime);
    fprintf('   Slip factor: \t\t%4.5e\n', sum.etac);
    fprintf('   Momentum compaction factor: \t %4.5e (%4.5e)\n',  sum.integrals(1)/sum.circumference, sum.compactionFactor);
    fprintf(SeparatorString);
    fprintf('   Assuming cavities Voltage: \t% 4.5f [kV]\n', v_cav/1e3);
    fprintf('                   Frequency: \t% 4.5f [MHz]\n', freq/1e6);
    fprintf('             Harmonic Number: \t% 4.0f\n', sum.harmon);
    fprintf('   Overvoltage factor: \t\t% 4.5f\n', sum.overvoltage);
    fprintf('   Synchronous Phase:  \t\t% 4.5f [rad] (%4.5f [deg])\n', sum.syncphase, sum.syncphase*180/pi);
    fprintf('   Linear Energy Acceptance:  \t% 4.3f %%\n', sum.energyacceptance*100);
    fprintf('   Synchrotron Tune:   \t\t% 4.5f (%4.5f kHz or %4.2f turns) \n', sum.synctune, sum.synctune/sum.revTime*1e-3, 1/sum.synctune);
    fprintf('   Bunch Length:       \t\t% 4.5f [mm] (%4.5f ps)\n', sum.bunchlength*1e3, sum.bunchlength/3e8*1e12);
    fprintf(SeparatorString);
    fprintf('   Injection:\n');
    fprintf('   H: beta = %06.3f [m] alpha = %+04.1e eta = %+04.3f [m] eta'' = %+04.1e \n', bx(1), ax(1), etax(1), etaprimex(1));
    fprintf('   V: beta = %06.3f [m] alpha = %+04.1e eta = %+04.3f [m] eta'' = %+04.1e \n', by(1), ay(1), etay(1), etaprimey(1));
    fprintf('   ********** End of Summary for ''%s'' **********\n', GLOBVAL.LatticeFile);
    fprintf('\n');
end

if nargout > 0
    varargout{1} = sum;
end

function [dI1,dI2,dI3,dI4,dI5,curHavg, Dxavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1,th1,th2)
%[dI1,dI2,dI3,dI4,dI5,curHavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1)
%calculate the contribution to the radiation integrals of a dipole.
%  INPUTS
%  rho, theta, radius and angle of the dipole
%  a0, b0, are horizontal alpha and beta at the entrance of the dipole,
%  D0, D0p are dispersion at the entrance of the dipole
%  K1, the gradient parameter in AT's convention, i.e., positive for
%  horizontal focusing, K1=0 by default
%  th1, th2, the entrance and exit angle, respectively, th1=th2= 0 [theta/2] by
%  default.
%

% If not combine dipole
if nargin == 6
    K1=0;
end

%
if nargin<9
    th1 = 0; %theta/2.0;
    th2 = 0; %theta/2.0;
end

% Edge focusing
M21 = tan(th1)/rho;
D0p = M21*D0+D0p;
a0 = -M21*b0+a0;

% split the dipole in N pieces
N = 100;
th = (0:N)/N*theta;

% Compute Twiss parameters inside dipole
for ii=1:length(th)
    [Dx(ii), Dxp(ii)] = calcdisp(rho, th(ii), D0, D0p, K1);
    [ax, bx] = calctwiss(rho, th(ii), a0, b0, K1);
    curHavg1(ii) = (Dx(ii)^2+(ax*Dx(ii)+bx*Dxp(ii))^2)/bx;
end

% Edge focusing
M21 = tan(th2)/rho;
Dxp(end) =  M21*Dx(end)+Dxp(end);
ax  = -M21*bx+ax;
curHavg1(end) = (Dx(end)^2+(ax*Dx(end)+bx*Dxp(end))^2)/bx;

% Average data
curHavg = ((curHavg1(1)+curHavg1(end))/2.0 + sum(curHavg1(2:end-1)))/(length(th)-1);
Dxavg   = ((Dx(1)+Dx(end))/2.0 + sum(Dx(2:end-1)))/(length(th)-1);

dI1 = ((Dx(1) + Dx(end))/2.0 + sum(Dx(2:end-1)))*theta/N;
dI2 = abs(theta/rho);
dI3 = abs(theta/rho^2);
dI4 = (1/rho^2 + 2*K1)*dI1  - (Dx(1)/rho^2*tan(th1) + Dx(end)/rho^2*tan(th2));
dI5 = curHavg*abs(theta/rho^2);

function [Dx, Dxp] = calcdisp(rho, theta, D0, D0p, K1)
%calcdisp - calculate dispersion function inside a combined-function dipole
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. D0 - Horizontal dispersion function at the entrance
%  4. D0p - DErivative of  Horizontal dispersion function at the entrance
%  5. K1 - Focusing
%
% Transfert matrix of A wedge dipole p58 Handbook af Accelerator Physics
s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Dx  =  D0*cos(sqK*s) + D0p/sqK*sin(sqK*s)+(1-cos(sqK*s))/rho/sqK^2;
    Dxp = -D0*sqK*sin(sqK*s)+D0p*cos(sqK*s)+sin(sqK*s)/rho/sqK;
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Dx =  D0*cosh(sqK*s) + D0p/sqK*sinh(sqK*s)+(-1+cosh(sqK*s))/rho/sqK^2;
    Dxp = D0*sqK*sinh(sqK*s)+D0p*cosh(sqK*s)+sinh(sqK*s)/rho/sqK;
    
end

function [ax, bx] = calctwiss(rho, theta, a0, b0, K1)
% calctwiss calculate twiss function inside a combined-function dipole manget
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. a0 - Horizontal alpha function at the entrance
%  4. b0 - Horizontal beta function at the entrance
%  5. K1 - Focusing
%
%  [beta ] = [  Mx11^2        -2*MX11*Mx12         Mx12^2   ] [beta0 ]
%  [alpha] = [ -Mx11*Mx21 Mx11*Mx22 + Mx11*Mx21   -Mx12*Mx22] [alpha0]
%  [gamma] = [  Mx21^2        -2*MX21*Mx22         Mx22^2   ] [gamma0]

Mx = calcMx(rho, K1,theta);
g0 = (1+a0^2)/b0;
twx2 = [Mx(1,1)^2, -2*Mx(1,1)*Mx(1,2), Mx(1,2)^2;
    -Mx(1,1)*Mx(2,1), Mx(1,1)*Mx(2,2)+Mx(1,2)*Mx(2,1),-Mx(1,2)*Mx(2,2);
    Mx(2,1)^2, -2*Mx(2,1)*Mx(2,2),Mx(2,2)^2] * [b0, a0, g0]';
ax = twx2(2);
bx = twx2(1);

function Mx = calcMx(rho,K1,theta)
% calcMx calculate transfer matrice of a combined-function dipole manget

s = rho*theta;
if K1>-1/rho^2 %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Mx = [cos(sqK*s), sin(sqK*s)/sqK; -sqK*sin(sqK*s), cos(sqK*s)];
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Mx = [cosh(sqK*s), sinh(sqK*s)/sqK; sqK*sinh(sqK*s), cosh(sqK*s)];
end
