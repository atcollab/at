% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/response matrix')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

% load lattice
s28d=load('/machfs/liuzzo/EBS/S28D/LATTICE/AT/S28Dmerged_PA.mat');

ring=s28d.LOW_EMIT_RING_INJ;
[l,t,c]=atlinopt(ring,0,1);
r0=ring;

% mark quadrupoles to use for tune matching
indqf1=find(atgetcells(ring,'FamName','QF1\w*'));
ring=atsetfieldvalues(ring,indqf1,'ForTuneF',1);                
indqd2=find(atgetcells(ring,'FamName','QD2\w*'));
ring=atsetfieldvalues(ring,indqd2,'ForTuneD',1);                

% set errors, large, AT does not find a closed orbit
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=5e-6*randn(size(ind));
dy=5e-6*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);

%% test tune matching
rerr=fittunedelta2fam(rerr,r0);