% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

% load lattice
load ESRFLattice.mat

%% get RM
speclab='DFSESRF';

modelrmfile=fullfile(pwd,['RMmodel' speclab '.mat']);%

if ~exist([modelrmfile],'file')
    
    ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        indSCor,...
        indQCor,...
        [],...
        [0 0 0 0 0 0]',...
        [1 2 3 7 8 9]); % orbit, dpp and dispersion rm to steerers
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end

% set errors
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=5e-6*randn(size(ind));
dy=5e-6*randn(size(ind));
dr=5e-6*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);
rerr=atsettilt(rerr,ind,dr);

%% apply correction
inCOD=[0 0 0 0 0 0]';
[l,~,~]=atlinopt(ring,0,indBPM);
refdispersion=zeros(2,length(indBPM));
refdispersion(1,:)=arrayfun(@(a)a.Dispersion(1),l);
refdispersion(2,:)=arrayfun(@(a)a.Dispersion(3),l);

% steerers limited, increase eigenvectors number
[rcor,inCOD,hs,vs]=atdispersionfreesteering(...
    rerr,...
    indBPM,...
    indHCor,...
    indVCor,...
    inCOD,...
    [... several correction iterations with different number of eigenvector
    [20 20];...
    [40 40];...
    [60 60];...
    [80 80];...
    [97 96];...
    [97 96];...
    [97 96]...
    ],...
    [true false],...
    1.0,...
    0.9,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    refdispersion,...
    [],...
    true);


o=findorbit6Err(rerr,indBPM,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. orbit');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. orbit');
saveas(gca,'OrbCor.fig');
 export_fig('OrbCor.jpg','-r300');

 
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
     delta=1e-4;
% get initial dispersion
o=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=finddispersion6Err(rcor,indBPM,indrfc,alpha,delta,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe-refdispersion(1,:),'.-');hold on; plot(sBPM,oxc-refdispersion(1,:),'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. disp');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. disp');
saveas(gca,'DispCor.fig');
 export_fig('DispCor.jpg','-r300');
 
 % steerers limited, increase eigenvectors number, no dispersion correction
[rcor,inCOD,hs,vs]=atdispersionfreesteering(...
    rerr,...
    indBPM,...
    indHCor,...
    indVCor,...
    inCOD,...
    [... several correction iterations with different number of eigenvector
    [20 20];...
    [40 40];...
    [60 60];...
    [80 80];...
    [97 96];...
    [97 96];...
    [97 96]...
    ],...
    [true false],...
    1.0,...
    0.0,...  dispersion weigth set to 0.0, no dispersion correction
    ModelRM,...
    zeros(2,length(indBPM)),...
    refdispersion,...
    [],...
    true);


o=findorbit6Err(rerr,indBPM,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. orbit');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. orbit');
saveas(gca,'OrbCorNoDisp.fig');
 export_fig('OrbCorNoDisp.jpg','-r300');

 
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
     delta=1e-4;
% get initial dispersion
o=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=finddispersion6Err(rcor,indBPM,indrfc,alpha,delta,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe-refdispersion(1,:),'.-');hold on; plot(sBPM,oxc-refdispersion(1,:),'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. disp');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. disp');
saveas(gca,'DispCorNoDisp.fig');
 export_fig('DispCorNoDisp.jpg','-r300');
 
 