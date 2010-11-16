function tauT = calc_TouschekPM(TD,dppPM,Trf,Ib,U0,coupling, sigE, emit_x)
%tauT = calc_TouschekPM(TD,dppPM,Trf,Ib,U0,coupling, sigE, emit_x)
%tauT = calc_TouschekPM(TD,dppPM,alpha,Ib,U0,coupling, sigE, emit_x)
%		Ib, mA, single bunch current
%		U0, MeV, one-turn energy loss
%       emit_x, nm-rad
%       coupling, average emit_y/emit_x
%       TD, lattice function structure (same as twissring output) 
%       dppPM, Nx2, positive/negative momentum aperture
%       Trf, a structure or a scalar, when a structure, it consists of [Frequency, HarmNumber, Energy, alpha,
%       Voltage], when a scalar it is alpha of the lattice (other
%       parameters are assumed default of SPEAR3).
%ex:    
%        load('dppAP_LAT_DW_withIDs_Mar20_07_normrf.mat')
%        [td, tune,chrom] = twissring(THERING,0,indextab, 'chrom', 1e-5);
%        dppPM = [deltap' deltam'];
%        tauT = calc_TouschekPM(td,dppPM,a1,100/280,1.04,0.064e-2, 0.001, 18)/3600; %hrs
%        disp(tauT)

e0 = 1.6e-19; %Coulomb
cspeed = 299792458; 
r0 = 2.82e-15; %m

U0=U0*1e6;
emit_x = emit_x*1.0e-9; %convert nm-rad to m-rad

%cavity related parameters
if isstruct(Trf)
    freq = Trf.Frequency;
    harm = Trf.HarmNumber;
    E0 = Trf.Energy;
    alpha = Trf.alpha;
    Vrf = Trf.Voltage;
    circ = Trf.circum;
else
    freq = 352.202e6; %Hz
    harm = 992;
    E0 = 6.04e9; %eV
    alpha = Trf;
    Vrf = 9e6;
    circ = 844.39;
    
end
gamma = E0/0.511e6;
N0 = 0.001/(freq/harm)/e0; %Number of particle per 1mA bunch. 


%bunch length
phi_s = asin(U0/Vrf);
nus = sqrt(harm*Vrf*alpha*cos(phi_s)/2/pi/E0);


sigZ = sigE/nus*harm*alpha/2/pi/freq*cspeed;

%rf bucket height
delta_max_rf = sqrt(2*U0/pi/alpha/harm/E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));
%---------------------------------

%beam size around the ring
%[td, tune,chrom] = twissring(THERING,0,1:length(THERING)+1, 'chrom', 1e-5);
td = TD; 
Dx = cat(2, td.Dispersion)';
betxy = cat(1, td.beta);
alfxy = cat(1, td.alpha);

spos = cat(1,td.SPos);

sigX = sqrt(betxy(:,1)*emit_x+Dx(:,1).^2*sigE^2);

sigY = sqrt(betxy(:,2)*emit_x*coupling);
sigXp = sqrt(emit_x*(1+alfxy(:,1).^2)./betxy(:,1)+Dx(:,2).^2*sigE^2);
%--------------------------------

curH = (Dx(:,1).^2 + (betxy(:,1).*Dx(:,2)+alfxy(:,1).*Dx(:,1)).^2)./betxy(:,1);

%delta_max_perp = hori_acceptance./sqrt(curH);
deltap = dppPM(:,1);
deltam = dppPM(:,2);

delta_maxp = min([deltap, ones(size(curH))*delta_max_rf]')';
delta_maxm = min([-deltam, ones(size(curH))*delta_max_rf]')';

xip = (delta_maxp/gamma.*betxy(:,1)./sigX).^2;
xim = (delta_maxm/gamma.*betxy(:,1)./sigX).^2;
Dvalp = funcD(xip);
Dvalm = funcD(xim);

ds = diff(spos);
n=1:length(ds);

avgfacp = sum(Dvalp(n)./sigX(n)./sigY(n)/sigZ./delta_maxp(n).^3.*ds)/circ;
avgfacm = sum(Dvalm(n)./sigX(n)./sigY(n)/sigZ./delta_maxm(n).^3.*ds)/circ;

lossrate = Ib*N0*r0^2*cspeed/8/gamma^2/pi*(avgfacp+avgfacm)/2.;
tauT = 1/lossrate;

if 0
   figure
   h = plot(spos, delta_maxp, spos, -delta_maxm); %delta_max_rf*ones(size(spos))); 
   set(gca,'fontsize', 16,'xlim',[0,240])
   xlabel('s (m)')
   ylabel('\delta_{max}')
   grid
   set(gca,'ylim',[-0.03,0.03]);
   
end

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

