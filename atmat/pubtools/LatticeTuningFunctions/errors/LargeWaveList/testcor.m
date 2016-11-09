clear all
close all

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/wave')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/random')
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/errordisplayfunctions');

% load lattice
load ../../ESRFLattice.mat
r0=ring;

% wave errors
ie=1;

wltouse=1:0.5:3;
amplx=0.6e-3;
amplY=0.6e-3;
amplpsi=0.6e-3;

W=findspos(r0,length(r0)+1)./wltouse;

A=amplx/length(W)*randn(size(W));
errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
errwavestruct(ie).type='x';
errwavestruct(ie).A=A(end:-1:1);
errwavestruct(ie).W=W;
ie=ie+1;

A=amplY/length(W)*randn(size(W));
errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
errwavestruct(ie).type='y';
errwavestruct(ie).A=A(end:-1:1);
errwavestruct(ie).W=W;
ie=ie+1;

A=amplpsi/length(W)*randn(size(W));
errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
errwavestruct(ie).type='psi';
errwavestruct(ie).A=A(end:-1:1);
errwavestruct(ie).W=W;
ie=ie+1;

magindex=arrayfun(@(a)a.indx,errwavestruct,'un',0);
type=arrayfun(@(a)a.type,errwavestruct,'un',0);
A=arrayfun(@(a)a.A,errwavestruct,'un',0);
W=arrayfun(@(a)a.W,errwavestruct,'un',0);

rerr=atsetwaveerrors(...
    r0,...
    magindex,...
    findcells(r0,'Class','Monitor'),...
    W,...
    A,...
    type);

    
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr,'comment',[],@pltmisalignments);
saveas(gca,'Wave.fig')
export_fig('Wave.jpg','-r300')

