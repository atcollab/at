% test errors and correction functions
close all
clear all

addpath('../../../LatticeTuningFunctions');
addpath('../../../LatticeTuningFunctions/correction/response matrix')
addpath('../../../LatticeTuningFunctions/correction/');
addpath('../../../LatticeTuningFunctions/errors/');

% load lattice

%% get RM
speclab='ESRF';

load('MeasurementESRF_SR_2018Feb20_resp1.mat');
r0 = atsetRFCavity(r0,8e6,0,992,0.0);
ring = r0;
rerr =atsetRFCavity(rfitope,8e6,0,992,0.0); %lattice fitted by CTRM application for later comparison

inCOD = [0 0 0 0 0 0]';
indBPM = find(atgetcells(r0,'Class','Monitor'))';
indhsel = find(atgetcells(r0,'Class','Sextupole'))';
indvsel = find(atgetcells(r0,'Class','Sextupole'))';

cor_in_sext =[ 2     4     6     9    11    13    16    18    20    23    25    27    30    32    34    37    39    41 44    46    48    51    53    55    58    60    62    65    67    69    72    74    76    79    81    83 86    88    90    93    95    97   100   102   104   107   109   111   114   116   118   121   123   125 128   130   132   135   137   139   142   144   146   149   151   153   156   158   160   163   165   167 170   172   174   177   179   181   184   186   188   191   193   195   198   200   202   205   207   209 212   214   216   219   221   223  ];
   
% smaller set of correctors
indHCor = indhsel(cor_in_sext(1:6:end)); 
indVCor = indHCor;

% mark fit locations
iq = find(atgetcells(r0,'Class','Quadrupole'));
ib = find(atgetcells(r0,'Class','Bend'));
% ind = [iq(1:40:end); ib(1:40:end)];
ind = [iq(1:1:end); ib([1:8:end,4:8:end,5:8:end,8:8:end])];
% ind = [iq(1:1:end); ib(1:end)];
r0 = atsetfieldvalues(r0,ind,'FitElement',1);

% measured RM structure
f0=992*PhysConstant.speed_of_light_in_vacuum.value/findspos(r0,length(r0)+1);

fresph=f0*oh;   % <---- qemcheckresp [m/Hz] * Hz
frespv=f0*ov;

MeasRMH = [resph(:);respv(:);oh(:);0.44;0.39];
MeasRMV = [resph2v(:);respv2h(:);ov(:)];

% structrue for measured RM
%  
delta=1e-4;
alpha=mcf(r0);
indrfc=find(atgetcells(r0,'Frequency'));

rmfunctnorm =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'normdisp',[0 0],txt); 
rmfunctskew =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'skewdisp',0,txt);
rmfunctdh =@(r,ib,~,~,~,~)getdisph6D(r,ib,indrfc,alpha,delta,inCOD);
rmfunctdv =@(r,ib,~,~,~,~)getdispv6D(r,ib,indrfc,alpha,delta,inCOD);

measuredRM.rmfunctnorm = rmfunctnorm;
measuredRM.rmfunctskew = rmfunctskew;
measuredRM.rmfunctdh = rmfunctdh;
measuredRM.rmfunctdv = rmfunctdv;
measuredRM.RMH = MeasRMH;%RMH;
measuredRM.RMV = MeasRMV;%RMV;
measuredRM.DH = oh;%DH;
measuredRM.DV = ov;%DV;

% expected perfect RM
RMH=rmfunctnorm(r0,indBPM,indHCor,indVCor,[],'computing hor. cor. RM');
RMV=rmfunctskew(r0,indBPM,indHCor,indVCor,[],'computing ver. cor. RM');
DH=rmfunctdh(r0,indBPM,[],[],[],'');
DV=rmfunctdv(r0,indBPM,[],[],[],'');
   
figure; plot(RMH(:)); hold on; plot(measuredRM.RMH(:)); legend('theory','measurement')
figure; plot(DH(:)); hold on; plot(measuredRM.DH(:)); legend('theory','measurement')
      
    
% number of eigenvectors for fit [quad, dip, skew]
neig = [100 64 100];

nsubsets = 4;

[...
    rfit,...        % fitted lattice
    Kqn,...         % quad gradient errors
    Kqs,...         % skew gradient errors
    Kdh,...         % bending angle errors
    Kdv,...         % dipole rotation errors
    indquad,...     % index of Kqn and Kqs
    inddip...       % index of Kdh Kdv
    ]=atFitResponseMatrixAndDispersion(...
    r0,...          1) lattice with errors to model (IGNORED)
    r0,...          2) lattice without errors
    'indBPM',indBPM,...     3) bpm indexes for rm fit
    'indHCor',indHCor,...     4) h correctors for rm fit
    'indVCor',indVCor,...     5) v correctors for rm fit
    'neig',neig,...           6) # eigenvectors for fit [quad, dip, skew]
    'nsubsets',nsubsets,...   7) # subsets for fit [quad, skew] errors=<fit(subsets)>
    'modecalc','Analytic',...
    'speclab',speclab,...      8) label to distinguish rm files and avoid recomputation if already existing
    'measuredRM',measuredRM);


%% plot fit result

alpha=mcf(r0);
indrfc=find(atgetcells(r0,'Frequency'));
delta =1e-4;
inCOD = zeros(6,1);

% reference lattice
[l,t,ch]=atlinopt(r0,0,indBPM);
bx0=arrayfun(@(a)a.beta(1),l);
by0=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(r0,indBPM,indrfc,alpha,delta,inCOD);
dx0=d(1,:);
dy0=d(3,:);

% lattice with errors
[l,t,ch]=atlinopt(rerr,0,indBPM);
bxe=arrayfun(@(a)a.beta(1),l);
bye=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
dxe=d(1,:);
dye=d(3,:);


[l,t,ch]=atlinopt(rfit,0,indBPM);
bxf=arrayfun(@(a)a.beta(1),l);
byf=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(rfit,indBPM,indrfc,alpha,delta,inCOD);
dxf=d(1,:);
dyf=d(3,:);

%
figure; 
plot((bxe - bx0)./bx0); hold on;
plot((bxf - bx0)./bx0); hold on;
xlabel('BPM #')
ylabel('\Delta\beta_x/\beta_{x,0}')
legend('errors','fit');
figure; 
plot((bye - by0)./by0); hold on;
plot((byf - by0)./by0); hold on;
xlabel('BPM #')
ylabel('\Delta\beta_y/\beta_{y,0}')
legend('errors','fit');

figure; 
plot((dxe - dx0)); hold on;
plot((dxf - dx0)); hold on;
xlabel('BPM #')
ylabel('\Delta\eta_x')
legend('errors','fit');
figure; 
plot((dye - dy0)); hold on;
plot((dyf - dy0)); hold on;
xlabel('BPM #')
ylabel('\Delta\eta_y')
legend('errors','fit');
