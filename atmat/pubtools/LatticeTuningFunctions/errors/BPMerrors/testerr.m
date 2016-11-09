clear all
close all

% load lattice
load ../../ESRFLattice.mat

% get indexes
indm=indBPM';
indq=find(atgetcells(ring,'Class','Quadrupole'));
%ring=atsetfieldvalues(ring,indm,'PassMethod','MonitorPass');

% define quadrupole alignment and rotation errors
dx=1e-6*randn(size(indq)); % random errors of 1um
dy=1e-6*randn(size(indq)); % random errors of 1um
dt=1e-6*randn(size(indq)); % random errors of 1urad

% define bpm offset and rotation errors
dox=1e-4*randn(size(indm)); % random misalignment errors at BPM of 100um
doy=1e-4*randn(size(indm)); % random misalignment errors at BPM of 100um
ox=1e-5*randn(size(indm)); % random offset errors of 10um
oy=1e-5*randn(size(indm)); 
gx=1e-3*randn(size(indm)); % random gain errors of 0.1%
gy=1e-3*randn(size(indm));  
rx=1e-6; % reading error sigma of 1um (can also be a vector)
ry=1e-6; 
rot=1e-5*randn(size(indm)); % random rotation errors of 10urad

% set errors
ringerr=ring;
%ringerr=atsetshift(ringerr,indq,dx,dy);
%ringerr=atsettilt(ringerr,indq,dt);
ringerr=atsetshift(ringerr,indm,dox,doy);
ringerr=atsetbpmerr(ringerr,indm,ox,oy,gx,gy,rx,ry,rot);

% plots
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
s=findspos(ringerr,indm);
o=findorbit4(ringerr,0,indm);
plot(s,o(1,:)'*1e6,'k');
hold on;
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
legend('orbit','bpm reading');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
xlabel('s [m]');ylabel('x [\mum]')
saveas(gca,'OrbitBPMAllErrX.fig')
export_fig('OrbitBPMAllErrX.jpg','-r300')

% plots
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
s=findspos(ringerr,indm);
o=findorbit4(ringerr,0,indm);
plot(s,o(1,:)'*1e6,'k');
hold on;
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
legend('orbit','bpm reading');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
oe=findorbit4Err(ringerr,0,indm);
plot(s,oe(1,:)'*1e6,'rx');
xlabel('s [m]');ylabel('y [\mum]')
saveas(gca,'OrbitBPMAllErrY.fig')
export_fig('OrbitBPMAllErrY.jpg','-r300')

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ringerr,[0,100],'comment',[],@plClosedOrbitBPM);
saveas(gca,'OrbitBPMAllErratplot.fig')
export_fig('OrbitBPMAllErratplot.jpg','-r300')

return
