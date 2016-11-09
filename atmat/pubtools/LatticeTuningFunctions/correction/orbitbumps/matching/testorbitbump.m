% test matching orbit bump
% test errors and correction functions
close all
clear all

addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/'));
addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/'));

% load lattice
s28d=load('/machfs/liuzzo/EBS/S28D/LATTICE/AT/S28Dmerged_PA.mat');

ring=s28d.LOW_EMIT_RING_INJ;
[l,t,c]=atlinopt(ring,0,1);
r0=ring; % lattice without errors


% get indexes
indBPM=find(atgetcells(ring,'Class','Monitor'))';
indHCor=find(atgetcells(ring,'iscorH','H'))';
indVCor=find(atgetcells(ring,'iscorV','V'))';
indSCor=find(atgetcells(ring,'iscorS','S'))';
indQCor=find(atgetcells(ring,'Class','Quadrupole'))';

inCOD=zeros(6,1);


bumph=1e-4;
bumpv=1e-4;
indBPMbump =indBPM([10]);
indHCorbump=indHCor([8 9 10]);
indVCorbump=indVCor([8 9 10]);
doplot=1;

[rbump,hs,vs]=BumpAtBPM(...
    ring,...
    inCOD,...
    bumph,...
    bumpv,...
    indBPMbump,...
    indHCorbump,...
    indVCorbump);

s=findspos(rbump,[indBPMbump,indHCorbump]);
figure;
atplot(rbump,[min(s)-2 max(s)+2],'comment',[],@plClosedOrbitOnly);
text(s(1)+0.1,bumph*1.01,{['hs:' num2str(hs')] ['vs:' num2str(vs') ]})
hold on; plot([s(1) s(1)],[0 bumph],'k:');
hold on; plot([s(1) s(1)],[0 bumpv],'r:');

%%
[rbump,hs,vs]=BumpAtBPM4D(...
    ring,...
    inCOD,...
    bumph,...
    bumpv,...
    indBPMbump,...
    indHCorbump,...
    indVCorbump);

s=findspos(rbump,[indBPMbump,indHCorbump]);
figure;
atplot(rbump,[min(s)-2 max(s)+2],'comment',[],@plClosedOrbitOnly);
text(s(1)+0.1,bumph*1.01,{['hs:' num2str(hs')] ['vs:' num2str(vs') ]})
hold on; plot([s(1) s(1)],[0 bumph],'k:');
hold on; plot([s(1) s(1)],[0 bumpv],'r:');
