
clear all
close all

% load lattice
load ../../ESRFLattice.mat

% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));

% set multipole errors
bn=[0 0 0 0 0 1.10934 0 0 0 -5.18658]*1e-4;

an=[
    0
    0
    4.803458371
    1.910276957
    1.055734675
    0.588073151
    0.312742308
    0.175288289
    0.101114708
    0.064747269]'*1e-4;

[rerr,PolB,PolA]=AssignFieldErr(ring,indq,2,7*1e-3,bn,an);


