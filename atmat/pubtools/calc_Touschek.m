function tauT = calc_Touschek(THERING,Ib,varargin)
%tauT = calc_Touschek(THERING, Ib)
%tauT = calc_Touschek(THERING, Ib,hori_acceptance)
%tauT = calc_Touschek(THERING, Ib,hori_acceptance,U0)
%tauT = calc_Touschek(THERING, Ib,hori_acceptance,U0,coupling)
%tauT = calc_Touschek(THERING, Ib,hori_acceptance,U0,coupling,sigE,emit_x)
%   hori_acceptance = Min(X/sqrt(beta)) around the ring
%		Ib, mA, single bunch current
%		Nb, number of bunches
%		U0, MeV, one-turn energy loss
%       emit_x, nm-rad
%ex: 
%tauT = calc_Touschek(THERING, 100/280, 0.015/sqrt(10.37),1.04,0.064e-2,0.001, 18)/3600
%
alpha=1.7e-4;

hori_acceptance = Inf;
if nargin>=3
	hori_acceptance = varargin{1};
end

if nargin<7
   atsum = atsummary;
end


if nargin>=4
	U0 = varargin{2}*1e6; %eV
else
    U0 = atsum.radiation*1e9; %eV
%U0 = 1.04e6; %eV
end
coupling = 0.05*1e-2; %by default
if nargin>=5
	coupling = varargin{3};
end

if nargin>=7
   sigE = varargin{4};
   emit_x = varargin{5}*1e-9; %m-rad
else
    sigE = atsum.naturalEnergySpread;  %sigma_delta
    emit_x = atsum.naturalEmittance;

end

e0 = PhysConstant.elementary_charge.value; %Coulomb
cspeed = PhysConstant.speed_of_light_in_vacuum.value; 
r0 = PhysConstant.classical_electron_radius.value; %m


%cavity related parameters
cava = findcells(THERING,'PassMethod','ThinCavityPass');
cavb = findcells(THERING,'PassMethod','CavityPass');
CAVINDEX = sort([cava,cavb]); %ati1.RF;
if isempty(CAVINDEX)
	error('cavity not defined')
end
freq = THERING{CAVINDEX(1)}.Frequency;
harm = THERING{CAVINDEX(1)}.HarmNumber;
E0 = THERING{CAVINDEX(1)}.Energy;
gamma = THERING{CAVINDEX(1)}.Energy/PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;

Vrf = 0;
for ii=1:length(CAVINDEX)
	Vrf = Vrf + THERING{CAVINDEX(ii)}.Voltage;
end
%Vrf = 3.2e6;
Vrf = 1.6e6;
%[alpha,a2] = findmcf3(THERING);

%bunch length
phi_s = asin(U0/Vrf);
nus = sqrt(harm*Vrf*alpha*cos(phi_s)/2/pi/E0);

sigZ = sigE/nus*harm*alpha/2/pi/freq*cspeed;

%rf bucket height
delta_max_rf = sqrt(2*U0/pi/alpha/harm/E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));
%---------------------------------

%beam size around the ring
[td, ~,~] = twissring(THERING,0,1:length(THERING)+1, 'chrom', 1e-5);
Dx = cat(2, td.Dispersion)';
betxy = cat(1, td.beta);
alfxy = cat(1, td.alpha);

spos = cat(1,td.SPos);
circ = spos(end);

sigX = sqrt(betxy(:,1)*emit_x+Dx(:,1).^2*sigE^2);

sigY = sqrt(betxy(:,2)*emit_x*coupling);
%sigXp = sqrt(emit_x*(1+alfxy(:,1).^2)./betxy(:,1)+Dx(:,2).^2*sigE^2);
%--------------------------------
curH = (Dx(:,1).^2 + (betxy(:,1).*Dx(:,2)+alfxy(:,1).*Dx(:,1)).^2)./betxy(:,1);

disp('delta_max_perp data:  ');
delta_max_perp = hori_acceptance./sqrt(curH);
disp('delta_max data:  ');
delta_max = min([delta_max_perp, ones(size(curH))*delta_max_rf]')';
disp('xi data:  ');
xi = (delta_max/gamma.*betxy(:,1)./sigX).^2;
Dval = funcD(xi);

N0 = 0.001/(freq/harm)/e0; %Number of particle per 1mA bunch. 

ds = diff(spos);
n=1:length(THERING);

avgfac = sum(Dval(n)./sigX(n)./sigY(n)/sigZ./delta_max(n).^3.*ds)/circ;
lossrate = Ib*N0*r0^2*cspeed/8/gamma^2/pi*avgfac;

tauT = 1/lossrate;

% if 0
%    figure
%    plot(spos, delta_max, spos, delta_max_rf*ones(size(spos))); 
%    %set(gca,'fontsize', 16,'xlim',[0,240])
%    set(gca,'fontsize', 16,'xlim',[0,120])
%    xlabel('s (m)')
%    ylabel('\delta_{max}')
%    grid
%    %set(gca,'ylim',[0,0.04]);
%    set(gca,'ylim',[0,0.15]);
%    
% end

function D=funcD(xi)
%a look-up table
DfunTable = [
%xi				Dfunc
0.000500	0.123802	
0.001000	0.153464	
0.001500	0.172578	
0.002000	0.186757	
0.002500	0.198008	
0.003000	0.207298	
0.003500	0.215179	
0.004000	0.221992	
0.004500	0.227968	
0.005000	0.233269	
0.005500	0.238015	
0.006000	0.242294	
0.006500	0.246176	
0.007000	0.249717	
0.007500	0.252961	
0.008000	0.255944	
0.008500	0.258697	
0.009000	0.261244	
0.009500	0.263607	
0.010000	0.265805	
0.010500	0.267852	
0.011000	0.269763	
0.011500	0.271549	
0.012000	0.273221	
0.012500	0.274788	
0.013000	0.276259	
0.013500	0.277640	
0.014000	0.278938	
0.014500	0.280159	
0.015000	0.281308	
0.015500	0.282391	
0.016000	0.283411	
0.016500	0.284372	
0.017000	0.285278	
0.017500	0.286132	
0.018000	0.286938	
0.018500	0.287698	
0.019000	0.288415	
0.019500	0.289090	
0.020000	0.289727	
0.020500	0.290327	
0.021000	0.290893	
0.021500	0.291425	
0.022000	0.291926	
0.022500	0.292397	
0.023000	0.292840	
0.023500	0.293256	
0.024000	0.293646	
0.024500	0.294011	
0.025000	0.294352	];
ximin = DfunTable(1,1);
ximax = DfunTable(end,1);
xi(find(xi<ximin)) = ximin;
xi(find(xi>ximax)) = ximax;
D = interp1(DfunTable(:,1), DfunTable(:,2), xi,'linear');

