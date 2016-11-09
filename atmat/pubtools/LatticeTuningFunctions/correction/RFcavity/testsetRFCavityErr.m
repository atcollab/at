% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/response matrix')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

% load lattice
load ESRFLattice.mat

%% 

inCOD=[0 0 0 0 0 0]';

speclab='RF';

modelrmfile=fullfile(pwd,['RMmodel' speclab '.mat']);%

if ~exist([modelrmfile],'file')
    
    ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        [],...
        [],...
        [],...
        [],...
        [],...
        inCOD,...
        [3]);
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end

% set errors, small, AT does find a closed orbit
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=1.0e-6*randn(size(ind));
dy=1.0e-6*randn(size(ind));
dr=1.0e-6*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);
rerr=atsettilt(rerr,ind,dr);


%% 
rfv=9.0e6; harm=992;

[rerr]=atsetRFCavity(rerr,rfv,1,harm,0);

[...
    rcor,...
    inCODcor,...
    fc....
    ]=atRFcorrection(...
    rerr,...
    indBPM,...
    inCOD,...
    [1 1 1 1 1],...
    1,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    true);


o=findorbit6Err(rerr,indBPM,inCOD);
oxe=o(5,:);
oye=o(6,:);

o=findorbit6Err(rcor,indBPM,inCODcor);
oxc=o(5,:);
oyc=o(6,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('energy deviation');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('time lag');
saveas(gca,'RFCor.fig');
export_fig('RFCor.jpg');

return