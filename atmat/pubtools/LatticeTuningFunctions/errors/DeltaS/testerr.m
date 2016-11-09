clear all
close all

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/random')
% load lattice
load ../../ESRFLattice.mat

% get indexes
indm=indBPM;
indq=find(atgetcells(ring,'Class','Quadrupole')); % girders are defined by GS and GE markers (start and end of girder)

rerr=atsetrandomerrors(...
    ring,...
    indq,...
    indm,...
    1,...
    5e-2,...
    2.5,...
    's');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ring)
hold on;
atplotsyn(gca,rerr);

saveas(gca,'DeltaSQuad.fig')
export_fig('DeltaSQuad.jpg','-r300')
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr)
%% get indexes
indm=find(atgetcells(ring,'Class','Monitor'));
indd=find(atgetcells(ring,'Class','Bend')); % girders are defined by GS and GE markers (start and end of girder)

rerr=atsetrandomerrors(...
    ring,...
    indd(1:4:end),...
    indm,...
    123456,...
    5e-2,...
    2.5,...
    's');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ring)
hold on;
atplotsyn(gca,rerr);

saveas(gca,'DeltaSDipZoom.fig')
export_fig('DeltaSDipZoom.jpg','-r300')

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr)
