%DEMOKNOB illustrates the use of MATLAB GUI controls with AT
%   to modify several parameters in THERING
%   In this case - quadrupole strength of QF and QD families 
%   This demo creates spear2 lattice
%   and 2 control slider: for QF and for QD quadrupole families

clear all
close all
spear2
qflist =  findcells(THERING,'FamName','QF');
qdlist =  findcells(THERING,'FamName','QD');
QFW =  0.02;
QDW = -0.02;

weights1 = num2cell(QFW*ones(size(qflist)));
kn1 = struct('Position',num2cell(qflist),'FieldName','K','M',1,'N',1,'Weight',weights1);
% When a value is changed, the PLOTBETA callback renews the beta-functions plot
atslider(kn1,'Demo: Adjust QF','plotbeta');
%atslider(kn1,'Demo: Adjust QF');

weights2 = num2cell(QDW*ones(size(qdlist)));
kn2 = struct('Position',num2cell(qdlist),'FieldName','K','M',1,'N',1,'Weight',weights2);
atslider(kn2,'Demo: Adjust QD','plotbeta');
%atslider(kn2,'Demo: Adjust QD');