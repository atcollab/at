% test errors and correction functions
close all
clear all

addpath('../../../LatticeTuningFunctions');
addpath('../../../LatticeTuningFunctions/correction/response matrix')
addpath('../../../LatticeTuningFunctions/correction/');
addpath('../../../LatticeTuningFunctions/errors/');

% load lattice
load ESRFLattice.mat
r0=ring;
%% get RM
speclab='ESRF';

% set errors
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=5e-6*randn(size(ind));
dy=5e-6*randn(size(ind));
dr=5e-6*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);
rerr=atsettilt(rerr,ind,dr);

inCOD = [0 0 0 0 0 0]';
indbsel = find(atgetcells(r0,'Class','Monitor'))';
indhsel = find(atgetcells(r0,'Class','Sextupole'))';
indvsel = find(atgetcells(r0,'Class','Sextupole'))';

% smaller set of correctors
indhsel = indhsel(1:14:end); indvsel = indhsel;

% mark fit locations
iq = find(atgetcells(r0,'Class','Quadrupole'));
ib = find(atgetcells(r0,'Class','Bend'));
% ind = [iq(1:40:end); ib(1:40:end)];
ind = [iq(1:1:end); ib(1:1:end)];
r0 = atsetfieldvalues(r0,ind,'FitElement',1);

% number of eigenvectors for fit [quad, dip, skew]
neig = [100 100 100];

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
    rerr,...        1) lattice with errors to model
    r0,...          2) lattice without errors
    inCOD,...       3) guess for initial coordinates
    indbsel,...     3) bpm indexes for rm fit
    indhsel,...     4) h correctors for rm fit
    indvsel,...     5) v correctors for rm fit
    neig,...        6) # eigenvectors for fit [quad, dip, skew]
    nsubsets,...    7) # subsets for fit [quad, skew] errors=<fit(subsets)>
    speclab...      8) label to distinguish rm files and avoid recomputation if already existing
    );


%% plot fit result
indBPM = indbsel;

alpha=mcf(r0);
indrfc=find(atgetcells(r0,'Frequency'));
delta =1e-4;
inCOD = zeros(6,1);

% reference lattice
[l,t,ch]=atlinopt(r0,0,indHCor);
bx0=arrayfun(@(a)a.beta(1),l);
by0=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(r0,indBPM,indrfc,alpha,delta,inCOD);
dx0=d(1,:);
dy0=d(3,:);

% lattice with errors
[l,t,ch]=atlinopt(rerr,0,indHCor);
bxe=arrayfun(@(a)a.beta(1),l);
bye=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
dxe=d(1,:);
dye=d(3,:);


[l,t,ch]=atlinopt(rfit,0,indHCor);
bxf=arrayfun(@(a)a.beta(1),l);
byf=arrayfun(@(a)a.beta(2),l);

d=finddispersion6Err(rfit,indBPM,indrfc,alpha,delta,inCOD);
dxf=d(1,:);
dyf=d(3,:);

%%
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



