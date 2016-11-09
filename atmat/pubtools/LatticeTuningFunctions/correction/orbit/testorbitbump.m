% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/response matrix')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

load ESRFLattice.mat

%% get RM
speclab='OrbitESRF';

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
        [1 2 3]);
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end

inCOD=[0 0 0 0 0 0]';

%% correct to reference orbit ( bump )

refx=zeros(size(indBPM));
refy=zeros(size(indBPM));

refx([7 8])=1e-4; % horizontal bump in position
refy([7 8])=[1e-4 -1e-4];  % vertical bump in angle

% get current orbit
o=findorbit6Err(ring,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

% correct restrictig response matrices
selbpm=[2  3    7 8   12 13]; % 9 and 12 missing in purpose
selcor=[2  3           4  5];

ModelRMbump.OrbHCor{1}=ModelRM.OrbHCor{1}(selbpm,selcor);
ModelRMbump.OrbVCor{3}=ModelRM.OrbVCor{3}(selbpm,selcor);
ModelRMbump.OrbHDPP=ModelRM.OrbHDPP;
ModelRMbump.OrbVDPP=ModelRM.OrbVDPP;
ModelRMbump.kval=ModelRM.kval;
ModelRMbump.delta=ModelRM.delta;

[rcor,inCOD]=atcorrectorbit(ring,...
    indBPM(selbpm),...
    indHCor(selcor),...
    indVCor(selcor),...
    inCOD,...
    repmat([4 4],3,1),...
    [false false],...
    1.0,...
    ModelRMbump,...
    [refx(selbpm);...
    refy(selbpm)],...
    [],...
    true);

o=findorbit6Err(rcor,indBPM,inCOD);
oxb=o(1,:);
oyb=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;
plot(sBPM,refx,':');hold on; plot(sBPM,refy,':');
plot(sBPM,oxb-oxc);hold on; plot(sBPM,oyb-oyc);
legend('ref hor.','ref ver.','hor.','ver.');
xlabel('BPM #');
ylabel('COD - COD_0')
xlim([0,100])
export_fig('CODbump.jpg')
%%
% recompute RM for subset of bpms and correctors.
[rcor,inCOD]=atcorrectorbit(ring,...
    indBPM(selbpm),...
    indHCor(selcor),...
    indVCor(selcor),...
    inCOD,...
    repmat([4 4],10,1),...
    [false true],... correctors average non zero and no frequency correction
    1.0,...
    [],...  <-- no RM, default is to compute it
    [refx(selbpm);...
    refy(selbpm)],...
    [],...
    true);


o=findorbit6Err(rcor,indBPM,inCOD);
oxb=o(1,:);
oyb=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;
plot(sBPM,refx,':');hold on; plot(sBPM,refy,':');
plot(sBPM,oxb-oxc);hold on; plot(sBPM,oyb-oyc);
legend('ref hor.','ref ver.','hor.','ver.');
xlabel('BPM #');
xlim([0,100])
ylabel('COD - COD_0')


%% set bump in lattice without errors
o=findorbit6Err(ring,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

% recompute RM for subset of bpms and correctors.
[rcor,inCOD]=atcorrectorbit(ring,...
    indBPM(selbpm),...
    indHCor(selcor),...
    indVCor(selcor),...
    inCOD,...
    repmat([4 4],5,1),...
    [false false],...
    1.0,...
    [],...  <-- no RM, default is to compute it
    [refx(selbpm);...
    refy(selbpm)],...
    [],...
    false);

o=findorbit6Err(rcor,indBPM,inCOD);
oxb=o(1,:);
oyb=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;
plot(sBPM,refx);hold on; plot(sBPM,refy);
plot(sBPM,oxb-oxc);hold on; plot(sBPM,oyb-oyc);
legend('ref hor.','ref ver.','hor.','ver.');
xlabel('BPM #');
ylabel('COD - COD_0')


