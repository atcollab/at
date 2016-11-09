% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');

addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/'));
addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/'));

% load lattice
load ESRFLattice.mat
speclab='ChainESRF';

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
        [],...
        [0 0 0 0 0 0]',...
        [1 2 3 4 5 6 7 8 9 10 11 12]); % all RM
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end


% mark quadrupoles to use for tune matching
indqf1=find(atgetcells(ring,'FamName','QF1\w*'));
ring=atsetfieldvalues(ring,indqf1,'ForTuneF',1);                
indqd2=find(atgetcells(ring,'FamName','QD2\w*'));
ring=atsetfieldvalues(ring,indqd2,'ForTuneD',1);                

inddq=find(atgetcells(ring,'FamName','DQ\w*'))';
inddl=find(atgetcells(ring,'FamName','DL\w*_3\w*'))';
ring=atsetfieldvalues(ring,[inddq inddl],'FitElement',1);     %mark as fitting point only some dipoles central ones.           

indDip=find(atgetcells(ring,'Class','Bend') & atgetcells(ring,'FitElement') )';

r0=ring;

% set errors, large, AT does not find a closed orbit
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=0.5e-4*randn(size(ind));
dy=0.5e-4*randn(size(ind));
dr=1.0e-4*randn(size(ind));
% truncate errors (quick wrong way, distribution not gaussian) 
dx(abs(dx)>2.5*0.5e-4)=0.5e-4;
dy(abs(dy)>2.5*0.5e-4)=0.5e-4;
dr(abs(dr)>2.5*1.0e-4)=1.0e-4;

rerr=atsetshift(ring,ind,dx,dy);
rerr=atsettilt(rerr,ind,dr);

%% correction chain

neigenvectors=[...
    200,... % n eig orbit H
    200,... % n eig orbit V
    200,... % skew quadrupole 
    250,... % normal quadrupole 
    250,... % fit normal quadrupole 
    100,... % fit dipole 
    250,... % fit skew quadrupole 
    ]; % number of eigenvectors 

diary('CorrChain7RF.txt');
cororder=[0 1 2 3 6 6 -1];
%  '(-1 ): RF cavity frequency and time lag tuning '...
%  '( 0 ): open trajectory (finds closed orbit) '...
%  '( 1 ): orbit '...
%  '( 2 ): tune '...
%  '( 3 ): chromaticity '...
%  '( 4 ): dispersion '...
%  '( 5 ): dispersion free steering '...
%  '( 6 ): rdt + dispersion correction '...

rcor=CorrectionChain(...
    rerr,...            %1  initial lattice
    r0,...              %2  model lattice
    indBPM,...          %3  bpm index
    indHCor,...         %4  h steerers index
    indVCor,...         %5  v steerers index
    indSCor,...  %6  skew quad index
    indQCor,...      %7  quadrupole correctors index
    neigenvectors,...            %8  number of eigen vectors [NeigorbitH, NeigorbitV, NeigQuadrdt, Neigdispv, Neigdisph,neig rdt corr, SkewQuadRDT]
    cororder,...       %9  correction order 1: orbit, 2: tune, 3: skewquad disp v 4: quad disp h 5: quad RDT 6: skew RDT
    ModelRM,...          %10 response matrices
    '',...          %11 response matrices
    true);

diary off


return
%%