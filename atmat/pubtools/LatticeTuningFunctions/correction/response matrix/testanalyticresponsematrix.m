% test errors and correction functions
close all
clear all

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/response matrix/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

% load lattice
s28d=load('/machfs/liuzzo/EBS/S28D/LATTICE/AT/S28Dmerged_PA.mat');

ring=s28d.LOW_EMIT_RING_INJ;
[l,t,c]=atlinopt(ring,0,1);
r0=ring;


speclab='FIT';

% get indexes
indBPM=find(atgetcells(ring,'Class','Monitor'))';
indHCor=find(atgetcells(ring,'iscorH','H'))';
indVCor=find(atgetcells(ring,'iscorV','V'))';
indSCor=find(atgetcells(ring,'iscorS','S'))';
indQCor=find(atgetcells(ring,'Class','Quadrupole'))';

% mark quadrupoles to use for tune matching
indqf1=find(atgetcells(ring,'FamName','QF1\w*'));
ring=atsetfieldvalues(ring,indqf1,'ForTuneF',1);                
indqd2=find(atgetcells(ring,'FamName','QD2\w*'));
ring=atsetfieldvalues(ring,indqd2,'ForTuneD',1);                

% %% lattice with errors

% set errors, small
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=0.5e-5*randn(size(ind));
dy=0.5e-5*randn(size(ind));
dr=5.0e-5*randn(size(ind));

rerr=ring;
% rerr=atsetshift(rerr,ind,dx,dy);
% rerr=atsettilt(rerr,ind,dr);
% 
% 
% % bpm resolution
% [l,t,ch]=atlinopt(r0,0,indBPM);
% Bx=arrayfun(@(a)a.beta(1),l);
% By=arrayfun(@(a)a.beta(2),l);
% bpmresx=(ones(size(indBPM))*1e-5./Bx)'; % better resolution for larger beam.
% bpmresy=(ones(size(indBPM))*1e-5./By)';
% % set BPM errors, from BPM fields. 
% rot=1e-4*randn(size(bpmresx));
% bpmoffx=0.0*bpmresx;
% bpmoffy=0.0*bpmresy;
% bpmgainx=bpmresx;
% bpmgainy=bpmresy;
% 
% rerr=atsetbpmerr(rerr,indBPM,bpmoffx,bpmoffy,bpmgainx,bpmgainy,bpmresx,bpmresy,rot);



% %% correct lattice if needed


inCOD=[0 0 0 0 0 0]';
% %% simulate RM measurement

indhsel=indHCor([1 2]);
indvsel=indVCor([1 2]);
indbsel=indBPM(:)';

delta=1e-4;
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
         
% "measure" RM, dispersion, tunes
rmfunctnorm =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'norm',[1e-2 1e-2],txt);
rmfunctskew =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'skew',1e-2,txt);
rmfunctdh =@(r,ib,~,~,~,~)getdisph6D(r,ib,indrfc,alpha,delta,inCOD);
rmfunctdv =@(r,ib,~,~,~,~)getdispv6D(r,ib,indrfc,alpha,delta,inCOD);

RMH=rmfunctnorm(rerr,indbsel,indhsel,indvsel,[],'computing hor. cor. RM');
RMV=rmfunctskew(rerr,indbsel,indhsel,indvsel,[],'computing ver. cor. RM');
DH=rmfunctnorm(rerr,indbsel,[],[],[],' ');
DV=rmfunctskew(rerr,indbsel,[],[],[],' ');


%% 

[rmh,rmv]=AnalyticSteerersResponse(rerr,indbsel,indhsel);

figure; plot(RMH); hold on; plot([rmh(:); rmv(:)]);%plot(RMH-[rmh(:); rmv(:)]);






