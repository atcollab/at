%FINDRESPMDEMO response matrix demo
% This script illustrates the use of AT function FINDRESPM

spear2
% The most common RM is corrector-to-BPM
% In this demonstration we will not use the actual correctors
% to keep the lattice simple. 

% We will use all focusing quadrupoles as correctors:
% In order to do this we need to use StrMPolesymplectic4 pass-method
% for them. This mehod looks at all terms of the polynomial
% expansion of transverse magnetic field. 
% (QuadLinearPass only looks at field 'K')
% PolynomB(1) gives horizontal kick
% PolynomA(1) gives a vertical kick

% Find indexes of elements that belong to QF Q1 Q2 Q3 families
% We will use them as corrector elements
QFI = findcells(THERING,'FamName','QF');
Q1I = findcells(THERING,'FamName','Q1');
Q2I = findcells(THERING,'FamName','Q2');
Q3I = findcells(THERING,'FamName','Q3');
CORRINDEX = sort([ QFI Q1I Q2I Q3I]);
% Install the new pass-method 'StrMPoleSymplectic4Pass'
THERING = setcellstruct(THERING,'PassMethod',CORRINDEX,'StrMPoleSymplectic4Pass');

% We will use etrance points of all bending magnets as observation points (BPMs)
BPMINDEX = findcells(THERING,'BendingAngle');

NBPM = length(BPMINDEX);
NCOR = length(CORRINDEX);

% Prepare input parameters for FINDRESPM that will tell it, which
% parameters to use as orbit perturbations
% See help for FINDRESPM

% Set the size of a parameter change for numeric differentiation
KICKSIZE = 1e-5;

RX = findrespm(THERING,BPMINDEX ,CORRINDEX, KICKSIZE, 'PolynomB',1,1);
RY = findrespm(THERING,BPMINDEX ,CORRINDEX, KICKSIZE, 'PolynomA',1,1);
% Build the response matrix 
% In the form
%
% | HH HV |
% | VH VV |
%
% HH - Horizontal BPM response to horizontal orbit kicks
% HV - Horizontal BPM response to vertical orbit kicks
% VH - vertical BPM response to horizontal orbit kicks
% VV - vertical BPM response to vertical orbit kicks
RespM_XY = [RX{1} RY{1}; RX{3} RY{3}];
figure(1);
mesh(RespM_XY);
colormap('copper');
xlabel('Corrector Number')
ylabel('BPM Number');
zlabel('Normalized Orbit Response');
title('Orbit Response Matrix - uncoupled lattice')

% Now we wish to introduce coupling:
QDI = findcells(THERING,'FamName','QD');
% Generate random rotations: 
QDTILTS =   1*(pi/180)*randn(1,length(QDI));
% Put random values in the ring
settilt(QDI,QDTILTS);

% Generate the new response matrix for the lattice with errors
RX = findrespm(THERING,BPMINDEX ,CORRINDEX, KICKSIZE, 'PolynomB',1,1);
RY = findrespm(THERING,BPMINDEX ,CORRINDEX, KICKSIZE, 'PolynomA',1,1);

RespM_XY_Coupled = [RX{1} RY{1}; RX{3} RY{3}];
figure(2);
mesh(RespM_XY_Coupled);
colormap('copper');
title('Orbit Response Matrix - coupled lattice')
xlabel('Corrector Number')
ylabel('BPM Number');
zlabel('Normalized Orbit Response');

