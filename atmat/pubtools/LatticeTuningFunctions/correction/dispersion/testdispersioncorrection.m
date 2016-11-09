% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

% load lattice
load ESRFLattice.mat

%% get RM
speclab='DispersionESRF';

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
        [9 10 11]);
    
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


[rcor,inCOD,hs,vs]=atcorrectdispersion(...
    rerr,...
    indBPM,...
    indQCor,...
    indSCor,...
    inCOD,...
    [floor(linspace(20,250,7));...
     floor(linspace(20,250,7))]',...
    [false true],...
    1.0,...
    ModelRM,...
    refdispersion,...
    [],...
    true);

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
ylabel('hor. disp - model');
subplot(2,1,2);
plot(sBPM,oye-refdispersion(2,:),'.-');hold on; plot(sBPM,oyc-refdispersion(2,:),'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. disp');
saveas(gca,'DispCor.fig');
export_fig('DispCor.jpg','-r300');
 
 return
