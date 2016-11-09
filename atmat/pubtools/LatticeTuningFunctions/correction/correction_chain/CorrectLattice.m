function rcor=CorrectLattice(...
    r0,...          1) lattice without error
    rerr,...        2) lattice to be corrected
    indBPM,...      3) monitor indexes in r0 and rerr
    indHCor,...     4) h steerers indexed in r0 and rerr
    indVCor,...     5) v steerers indexed in r0 and rerr
    indSCor,...     6) normal quad correctors indexed in r0 and rerr
    indQCor,...     7) skew quad correctors indexed in r0 and rerr
    ModelRM,...     8) response matrice structure, if [], rm are computed
    corparam,...    9) correction parameters, order, eigenvectors #
    speclab,...
    verbose)
%
% performs a loop of corrections as described in corparam
%
% rcor=correctlattice(...
%     r0,...          1) lattice without error
%     rerr,...        2) lattice to be corrected
%     indBPM,...      3) monitor indexes in r0 and rerr
%     indHCor,...     4) h steerers indexed in r0 and rerr
%     indVCor,...     5) v steerers indexed in r0 and rerr
%     indSCor,...     6) normal quad correctors indexed in r0 and rerr
%     indQCor,...     7) skew quad correctors indexed in r0 and rerr
%     ModelRM,...     8) response matrice structure, if [], rm are computed
%     corparam...     9) correction parameters, order, eigenvectors #
%     speclab,...    10) name for rm files default : ''
%     verbose...     11) if true print correction parameters at all steps
%     )
%
% corparam is a structure with fields
%
% corparam.cororder=[0,1,1,2,3,7,1,2,3,7];
%
%  '( 0): open trajectory (finds closed orbit) '...
%  '( 1): orbit '...
%  '( 2): tune '...
%  '( 3): chromaticity '...
%  '( 4): dispersion '...
%  '( 5): dispersion free steering '...
%  '( 6): rdt + dispersion correction '...
%  '( 7): fit errors model and correct model quad RDT + dispersion (6) '
%                
% corparam.neigeachcor=[...
%         200,...100,... % n eig orbit H
%         200,...100,... % n eig orbit V
%         100,...100,... % skew quadrupole vertical dispersion correction
%         100,...100,... % quadrupole horizontal dispersion correction
%         100,...100,... % quadrupole horizontal dispersion correction
%         100,...100,... % quadrupole horizontal dispersion correction
%         100,...100,... % quadrupole horizontal dispersion correction
%         ]; % number of eigenvectors orbit, quadrdt, disp v disp h skewRDT
%
%see also: CorrectionChain

if nargin<11 %
    verbose=false;
end
if nargin<10 %
    speclab='';
end
if nargin<9 %|| isempty(corparam)
    corparam.cororder=[0,1,1,2,3,7,1,2,3,7];
    
    corparam.neigenvectors=[...
        100,... % n eig orbit H
        100,... % n eig orbit V
        100,... % skew quadrupole vertical dispersion correction
        100,... % quadrupole horizontal dispersion correction
        100,... % n eig quad error fit
        100,... % n eig dipole error fit
        100,... % n eig skew error fit
        ]; % number of eigenvectors orbit, quadrdt, disp v disp h skewRDT
    
    
else
    cororder=corparam.cororder;
    neig=corparam.neigenvectors;
end

disp(['corr order  : ' num2str(cororder,'%d, ')]);
disp(['corr eig    : ' num2str(neig,'%d, ')]);

% correct
[rcor]=CorrectionChain(...
    rerr,...
    r0,...
    indBPM,...
    indHCor,...
    indVCor,...
    indSCor,...
    indQCor,...
    neig,...
    cororder,...
    ModelRM,...
    speclab,...
    verbose);
