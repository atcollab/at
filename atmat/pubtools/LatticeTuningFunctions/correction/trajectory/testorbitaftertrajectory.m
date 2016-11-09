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

%% get RM
speclab='OrbitAfterTrajectory';

% get indexes
indBPM=find(atgetcells(ring,'Class','Monitor'))';
indHCor=find(atgetcells(ring,'iscorH','H'))';
indVCor=find(atgetcells(ring,'iscorV','V'))';
indSCor=find(atgetcells(ring,'iscorS','S'))';
indQCor=find(atgetcells(ring,'Class','Quadrupole'))';

modelrmfile=fullfile(pwd,['RMmodel' speclab '.mat']);%

if ~exist([modelrmfile],'file')
    
    ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        indSCor,...
        indQCor,...
        [0 0 0 0 0 0]',...
        [1 2 3 4 5 6]);
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end

% set errors, large, AT does not find a closed orbit
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=1e-4*randn(size(ind));
dy=1e-4*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);


%%
corparam.neigenvectors=[...
    200,... % n eig orbit H
    200,... % n eig orbit V
    200,... % skew quadrupole 
    250,... % normal quadrupole 
    350,... % fit normal quadrupole 
    100,... % fit dipole 
    350,... % fit skew quadrupole 
    ]; % number of eigenvectors 

diary('CorrChainOrbitAfterTraj.txt');

corparam.cororder=[0 1];

rcor=CorrectLattice(...
    r0,... no errors
    rerr,... lattice to correct
    indBPM,...
    indHCor,...
    indVCor,...
    indSCor,...
    indQCor,...
    ModelRM,corparam,'');

diary off