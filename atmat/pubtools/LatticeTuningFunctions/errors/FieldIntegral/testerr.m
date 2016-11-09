clear all
close all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors')

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
    1e-1,...
    2.5,...
    'dpb2');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr,[0,100],'comment',[],@plPolynomBComp);
saveas(gca,'FieldIntegral.fig')
export_fig('FieldIntegral.jpg','-r300')

