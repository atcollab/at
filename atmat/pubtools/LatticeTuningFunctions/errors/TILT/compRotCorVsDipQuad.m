clear all
close all
addpath('/mntdirect/_users/liuzzo/atwork-code/physics/errorsettingfunctionsAT2')
% load lattice
load /machfs/liuzzo/ESRF/StorageRing/LATTICE/AT/ESRF_fortol.mat

% get indexes
indc=find(atgetcells(ring,'iscorH','H'));
indc=indc(45); % select a specific magnet for the test
ring0=ring;

%set correctors and dipole to 1e-5 rad horizontal kick and 1e-2 quad kick
kick=1e-5;
quad=0.01;
ringcor=ring;
ringcor=atsetfieldvalues(ringcor,indc,'PolynomB',{1,1},kick);
ringcor=atsetfieldvalues(ringcor,indc,'PolynomB',{1,2},quad);

br=getBrho(ring);
L=ring{indc}.Length;
ringdip=ring;
ringdip=atsetfieldvalues(ringdip,indc,'BendingAngle',{1,1},kick*L);
ringdip=atsetfieldvalues(ringdip,indc,'EntranceAngle',{1,1},0);
ringdip=atsetfieldvalues(ringdip,indc,'ExitAngle',{1,1},0);
ringdip=atsetfieldvalues(ringdip,indc,'PassMethod','BndMPoleSymplectic4FrgFPass');
ringdip=atsetfieldvalues(ringdip,indc,'PolynomB',{1,2},quad);


% define s-axis rotation errors in correctors
dt=1e-4; % rotation error of 1mrad

% set errors
ringcortilt=atsettilt(ringcor,indc,dt); % corrector
ringdiptiltref=atsettilt(ringdip,indc,dt); % dipole tilted
ringdiptilt=atsettiltdipole(ringdip,indc,dt); % dipole tilted
%ringdiptilt=settilt_THERING_Dipole(ringdip,indc,dt); % dipole tilted

% plots
[twin,~,~]=atlinopt(ring0,0,1);
indall=1:length(ring)+1;
lcor=twissline(ringcor,0.0,twin,indall,'chrom'); 
ldip=twissline(ringdip,0.0,twin,indall,'chrom'); 
lcortilt=twissline(ringcortilt,0.0,twin,indall,'chrom'); 
ldiptilt=twissline(ringdiptilt,0.0,twin,indall,'chrom'); 
ldiptiltref=twissline(ringdiptiltref,0.0,twin,indall,'chrom'); 
  
oxcor=arrayfun(@(a)a.ClosedOrbit(1),lcor,'un',1);
oycor=arrayfun(@(a)a.ClosedOrbit(3),lcor,'un',1);
dxcor=arrayfun(@(a)a.Dispersion(1),lcor,'un',1);
dycor=arrayfun(@(a)a.Dispersion(3),lcor,'un',1);

oxdip=arrayfun(@(a)a.ClosedOrbit(1),ldip,'un',1);
oydip=arrayfun(@(a)a.ClosedOrbit(3),ldip,'un',1);
dxdip=arrayfun(@(a)a.Dispersion(1),ldip,'un',1);
dydip=arrayfun(@(a)a.Dispersion(3),ldip,'un',1);
  
oxcortilt=arrayfun(@(a)a.ClosedOrbit(1),lcortilt,'un',1);
oycortilt=arrayfun(@(a)a.ClosedOrbit(3),lcortilt,'un',1);
dxcortilt=arrayfun(@(a)a.Dispersion(1),lcortilt,'un',1);
dycortilt=arrayfun(@(a)a.Dispersion(3),lcortilt,'un',1);

oxdiptilt=arrayfun(@(a)a.ClosedOrbit(1),ldiptilt,'un',1);
oydiptilt=arrayfun(@(a)a.ClosedOrbit(3),ldiptilt,'un',1);
dxdiptilt=arrayfun(@(a)a.Dispersion(1),ldiptilt,'un',1);
dydiptilt=arrayfun(@(a)a.Dispersion(3),ldiptilt,'un',1);

oxdiptiltref=arrayfun(@(a)a.ClosedOrbit(1),ldiptiltref,'un',1);
oydiptiltref=arrayfun(@(a)a.ClosedOrbit(3),ldiptiltref,'un',1);
dxdiptiltref=arrayfun(@(a)a.Dispersion(1),ldiptiltref,'un',1);
dydiptiltref=arrayfun(@(a)a.Dispersion(3),ldiptiltref,'un',1);


figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispCorQuad')
s=findspos(ring,indall);
yyaxis left
plot(s,oxcor,'b');hold on;
plot(s,oycor,'r');
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxcor,'c');hold on;
plot(s,dycor,'m');
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispCorQuad.fig')
export_fig('OrbitDispCorQuad.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuad')
s=findspos(ring,indall);
yyaxis left
plot(s,oxdip,'b');hold on;
plot(s,oydip,'r');
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxdip,'c');hold on;
plot(s,dydip,'m');
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuad.fig')
export_fig('OrbitDispDipQuad.jpg','-r300')


figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispCorQuadTilt')
s=findspos(ring,indall);
yyaxis left
plot(s,oxcortilt,'b');hold on;
plot(s,oycortilt,'r');
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxcortilt,'c');hold on;
plot(s,dycortilt,'m');
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispCorQuadTilt.fig')
export_fig('OrbitDispCorQuadTilt.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuadTilt')
s=findspos(ring,indall);
yyaxis left
plot(s,oxdiptilt,'b');hold on;
plot(s,oydiptilt,'r');
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxdiptilt,'c');hold on;
plot(s,dydiptilt,'m');
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuadTilt.fig')
export_fig('OrbitDispDipQuadTilt.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuadTiltRef')
s=findspos(ring,indall);
yyaxis left
plot(s,oxdiptiltref,'b');hold on;
plot(s,oydiptiltref,'r');
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxdiptiltref,'c');hold on;
plot(s,dydiptiltref,'m');
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuadTiltRef.fig')
export_fig('OrbitDispDipQuadTiltRef.jpg','-r300')


figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispCorQuadTiltVar')
s=findspos(ring,indall);
yyaxis left
plot(s,oxcortilt-oxcor,'b');hold on;
plot(s,oycortilt-oycor,'r');
ax=gca;
ax1ylim=ax.YLim;
xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,dxcortilt-dxcor,'c');hold on;
plot(s,dycortilt-dycor,'m');
ax=gca;
ax2ylim=ax.YLim;
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispCorQuadTiltVar.fig')
export_fig('OrbitDispCorQuadTiltVar.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuadTiltVar')
s=findspos(ring,indall);
yyaxis left
plot(s,oxdiptilt-oxdip,'b');hold on;
plot(s,oydiptilt-oydip,'r');
xlabel('s [m]');
ylabel('orbit [m]');
ax=gca;
ax.YLim=ax1ylim;

yyaxis right
plot(s,dxdiptilt-dxdip,'c');hold on;
plot(s,dydiptilt-dydip,'m');
ax=gca;
ax.YLim=ax2ylim;

xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuadTiltVar.fig')
export_fig('OrbitDispDipQuadTiltVar.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuadTiltVarRef')
s=findspos(ring,indall);
yyaxis left
plot(s,oxdiptiltref-oxdip,'b');hold on;
plot(s,oydiptiltref-oydip,'r');
xlabel('s [m]');
ylabel('orbit [m]');
ax=gca;
ax.YLim=ax1ylim;

yyaxis right
plot(s,dxdiptiltref-dxdip,'c');hold on;
plot(s,dydiptiltref-dydip,'m');
ax=gca;
ax.YLim=ax2ylim;

xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuadTiltVarRef.fig')
export_fig('OrbitDispDipQuadTiltVarRef.jpg','-r300')


figure('units','normalized','position',[0.3 0.3 0.45 0.35],'name','OrbitDispDipQuadVsTiltVar')
s=findspos(ring,indall);
yyaxis left
plot(s,(oxdiptilt-oxdip)-(oxcortilt-oxcor),'b');hold on;
plot(s,(oydiptilt-oydip)-(oycortilt-oycor),'r');
ax=gca;
ax.YLim=ax1ylim;

xlabel('s [m]');
ylabel('orbit [m]');
yyaxis right
plot(s,(dxdiptilt-dxdip)-(dxcortilt-dxcor),'c');hold on;
plot(s,(dydiptilt-dydip)-(dycortilt-dycor),'m');
ax=gca;
ax.YLim=ax2ylim;
xlabel('s [m]');
ylabel('dispersion [m]');
legend('x','y','\eta_x','\eta_y')
saveas(gca,'OrbitDispDipQuadVsTiltVar.fig')
export_fig('OrbitDispDipQuadVsTiltVar.jpg','-r300')


return



