clear all
close all

% load lattice
load ../../ESRFLattice.mat

% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));

% define s-axis rotation errors
dt=1e-6*randn(size(indq)); % random errors of 1um

% set errors
ringerr=atsettilt(ring,indq,dt);

% plots
figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ring,'comment',[],@plClosedOrbit)
saveas(gca,'OrbitNoErrTilt.fig')
export_fig('OrbitNoErrTilt.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],@plClosedOrbit)
saveas(gca,'OrbitWithErrTilt.fig')
export_fig('OrbitWithErrTilt.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
s=findspos(ring,indq);
plot(s,dt,'b.');
atplotsyn(gca,ringerr);
axis tight;
xlabel('s [m]');
ylabel('s-axis rotation [rad]');

saveas(gca,'SetErrTilt.fig')
export_fig('SetErrTilt.jpg','-r300')

