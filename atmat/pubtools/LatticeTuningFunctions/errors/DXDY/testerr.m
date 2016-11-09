clear all
close all

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/errordisplayfunctions');

% load lattice
load ../../ESRFLattice.mat

% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));

% define alignemnt errors
dx=1e-6*randn(size(indq)); % random errors of 1um
dy=1e-6*randn(size(indq));

% set errors
ringerr=atsetshift(ring,indq,dx,dy);

% plots
figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ring,'comment',[],@plClosedOrbit)
saveas(gca,'OrbitNoErr.fig')
export_fig('OrbitNoErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],@plClosedOrbit)
saveas(gca,'OrbitWithErr.fig')
export_fig('OrbitWithErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],@pltmisalignments)
saveas(gca,'SetErrDxDy.fig')
export_fig('SetErrDxDy.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,[0,52],'comment',[],@pltmisalignments)
saveas(gca,'SetErrDxDyZoom.fig')
export_fig('SetErrDxDyZoom.jpg','-r300')

