
close all
clear all
%% simple lattice test (uncomment to run)
% % file to be converted
% seqfilemadX='dbatest.seq';
% E0=0.51e9;
% MADX_2_AT(seqfilemadX,E0);
% 
% load([seqfilemadX(1:end-4) '_AT_LATTICE.mat']);
% % the variable RING is contained in the dbatest.seq lattice in MADX
% atplot(RING);

addpath(fullfile(pwd,'..'));
%% ESRF upgrade S10E test
seqfilemadX='low_emit_s10E_save.seq';
E0=6.039e9;

% execute conversion. file ..._AT_LATTICE.mat will be created.
% atfrommadx(seqfilemadX,E0,'anicename');
% load('anicename.mat')

atfrommadx(seqfilemadX,E0);
load([seqfilemadX(1:end-4) '_AT_LATTICE.mat']); 

atplot(ARC2);
