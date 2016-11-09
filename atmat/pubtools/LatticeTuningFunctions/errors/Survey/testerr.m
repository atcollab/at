clear all
close all

% load lattice
load /machfs/liuzzo/ESRF/StorageRing/LATTICE/AT/ESRF_fortol.mat



algedir='/machfs/liuzzo/EBS';
algeactfile=fullfile(algedir,'Actual_Position_Simu.xlsx');
rerr=SetESRFAlgeAlignmentError(ring,algeactfile,'',1);
    
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr,'comment',[],@pltmisalignments)

saveas(gca,'SurveyErrors.fig')
export_fig('SurveyErrors.jpg','-r300')
  
indm=find(atgetcells(rerr,'Class','Monitor'));

%% plots
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
s=findspos(rerr,indm);
o=findorbit4(rerr,0,indm);
plot(s,o(1,:)'*1e3,'k');
hold on;
oe=findorbit4Err(rerr,0,indm);
plot(s,oe(1,:)'*1e3,'rx');
legend('orbit','bpm reading');
oe=findorbit4Err(rerr,0,indm);
plot(s,oe(1,:)'*1e3,'rx');
oe=findorbit4Err(rerr,0,indm);
plot(s,oe(1,:)'*1e3,'rx');
oe=findorbit4Err(rerr,0,indm);
plot(s,oe(1,:)'*1e3,'rx');
oe=findorbit4Err(rerr,0,indm);
plot(s,oe(1,:)'*1e3,'rx');
xlabel('s [m]');ylabel('x [mm]')
saveas(gca,'OrbitBPMAllErrX.fig')
export_fig('OrbitBPMAllErrX.jpg','-r300')