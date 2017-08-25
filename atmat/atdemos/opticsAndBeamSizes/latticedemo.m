%LATTICEDEMO self-running tutorial
% demonstrates 
% 1. ELEMENT, family of ELEMENTS, sequence of ELEMENTS 
% 2. Lattice representation
% 3. Creating a lattice 
% 4. Creating and editing lattice files
clear all
clc
echo on
% An element in Accelerator Toolbox is a 1-by-1 MATLAB STRUCTURE
% The folowing code creates a structure D1 for a drift space
% and a structure QF for a quadrupole.

D1.FamName = 'DR01';
D1.Length  = 3;
D1.PassMethod = 'DriftPass';

QF.FamName = 'QF';
QF.Length = 1;
QF.K = 0.2;
QF.MaxOrder = 3;
QF.NumIntSteps = 1;
QF.PolynomA= [0 0 0];
QF.PolynomB= [0 0.2 0];
QF.R1= eye(6);
QF.R2= eye(6);
QF.T1= [0 0 0 0 0 0];
QF.T2= [0 0 0 0 0 0];
QF.PassMethod= 'QuadLinearPass';

pause % Press any key to continue 
clc
% Use WHOS, DISP or just type variable's name without closing semicolon 
% to print the element's info:

whos D1 QF                   
disp(D1)                      
QF                            

pause % Press any key to continue 
clc
% The next few lines will create another drift structure D2 from the exiting D1
% and modify the values of fields 'FamName' and 'Length'

D2 = D1;

D2.FamName = 'DR02';
D2.Length = 2;

disp(D2)
pause % Press any key to continue 
clc
% Create another quadrupole element structure QD from QF and modify
% the values of fields 'K' and 'PolynomB' to make it defocusing 
QD = QF;
QD.FamName = 'QD';
QD.K = -0.4;
% Field 'PolynomB is a vector with polynomial field expansion coefficients.
% The second element (quadrupole coefficient) must be consistent with field 'K' 
QD.PolynomB(2) = QD.K;

disp(QD)
pause % Press any key to continue 
clc
% We have declared four elements:
whos

% They are now independent from each other

% We are ultimately interested in sequences of elements
% to model storage ring lattices or single-pass beam transport lines.
% The next section will illustrate building of such sequences



pause % Press any key to continue 
clc
% Accelerator Toolbox represents sequences of elements as MATLAB cell arrays
% where individual cells are 1-by-1 structures containing element data
% The following commad creates a simple FODO cell by copying previously 
% created element structures for drifts and quadrupole magnets to a cell array FODOCELL:

FODOCELL = {QF D1 QD D2 QF};

whos FODOCELL
% LENGTH is useful to find the number of elements in a sequence

L = length(FODOCELL) 
pause % Press any key to continue;
clc
% Use {:} cell array syntax to print some or all elements
FODOCELL{1}
pause % FODOCELL{:} will print a long list of all elements. Press any key to continue
clc
FODOCELL{:}
pause % Press any key to continue;
clc
% Let's build a cell array THERING that represents a closed ring 
% with 10 periods of FODOCELL the same way we would build 
% any other array in MATLAB from the comman dline

THERING = [FODOCELL FODOCELL FODOCELL FODOCELL FODOCELL...
           FODOCELL FODOCELL FODOCELL FODOCELL FODOCELL]; 
        
whos THERING
pause % Press any key to continue;
clc
% The first element in THERING is 
THERING{1}

% To inspect or change the value of a specific field we can use MATLAB syntax
% for accessing cells in cell arrays and field in structures
oldK = THERING{1}.K

THERING{1}.K = 0.25;

newK = THERING{1}.K

pause % Press any key to continue;
clc
% Lattice THERING is a variable in MATLAB workspace.
% We can use it in accelerator physics functions and scripts
%
% For example: function FindM44 finds 4-by-4 transverse transfer matrix
M = findm44(THERING,0)
pause % Press any key to continue;
clc
% -----------------------------------------------------------------------
% SUMMARY
% 1. Individual elements are represented by 1-by-1  MATLAB structures

% 2. Element sequences (lattices) are represented by 1-dimensional 
%    MATLAB cell arrays of stuctures

% 3. MATLAB syntax for hanling structures and cell arrays applies. 
%    No special language is required to define a lattice. 

% --------------------------------------------------------------------------

pause % Press any key to continue;
echo off
clc
