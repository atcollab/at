clear all
close all
addpath('../../errors')
addpath('../../../../../machine_data')

% load lattice
ring=esrf;

% get indexes
indm=find(atgetcells(ring,'Class','Monitor'));
indq=find(atgetcells(ring,'Class','Quadrupole')); % girders are defined by GS and GE markers (start and end of girder)

rerr=atsetrandomerrors(...
    ring,...
    indq,...
    indm,...
    1,...
    1e-3,...
    2.5,...
    'dpb2');

figure('units','normalized','position',[0.1 0.4 0.65 0.45])
atplot(rerr,'comment',[],@plPolynomBComp);
figure('units','normalized','position',[0.1 0.4 0.65 0.45])
atplot(ring);
figure('units','normalized','position',[0.1 0.4 0.65 0.45])
atplot(rerr);
%saveas(gca,'FieldIntegral.fig')
%export_fig('FieldIntegral.jpg','-r300')

