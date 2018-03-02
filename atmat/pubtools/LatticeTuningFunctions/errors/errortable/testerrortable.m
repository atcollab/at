% test the possibility to use matlab tables for error handling
clear all
close all

% load lattice
addpath('/mntdirect/_users/liuzzo/ATcode/at/machine_data');
r = esrf;

% create table of errors
ErrTab = atcreateerrortable(r);

% enter some errors in the table

% hor misal
ErrTab.X = randn(size(ErrTab.X)).*1e-6 ;

% roll
ErrTab.Roll = randn(size(ErrTab.X)).*1e-6 ;

% quadrupole errors
mask = atgetcells(r,'Class','Quadrupole');

ErrTab.DK_K(mask,2) = randn(size(ErrTab.DK_K(mask,2))).*1e-4 ;

% Yaw
ErrTab.Yaw(mask) = randn(size(ErrTab.X(mask))).*1e-3 ;

% Pitch
ErrTab.Pitch(mask) = randn(size(ErrTab.X(mask))).*1e-3 ;

% BPM errors
mask = atgetcells(r,'Class','Monitor');
ErrTab.BPM_Offset(mask,:) = randn(length(ErrTab.X(mask)),2).*1e-4 ;
ErrTab.BPM_Reading(mask,:) = randn(length(ErrTab.X(mask)),2).*1e-6 ;

% systematic multipoles
mask = atgetcells(r,'Class','Sextupole');
ErrTab.b_n_systematic(mask,:) = repmat([0 0 0 0 1e-1 1e-2 1e-3 zeros(1,13)],length(find(mask)),1) ;

% random multipoles
mask = atgetcells(r,'Class','Sextupole');
ErrTab.a_n_random(mask,:) = repmat([0 0 0 0 1e-1 1e-2 1e-3 zeros(1,13)],length(find(mask)),1).*randn(length(find(mask)),20) ;

% apply them in the lattice
rerr = atseterrortable(r,ErrTab,'verbose',true);



