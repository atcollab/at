%LATTICEDEMO self-running tutorial
% demonstrates:
% 1. ELEMENT, family of ELEMENTS, sequence of ELEMENTS 
% 2. Lattice representation
% 3. Creating a lattice 
% 4. Creating and editing lattice files
clear
clc
echo on
% An element in Accelerator Toolbox is a 1-by-1 MATLAB STRUCTURE.
% Functions are provided to easuly create such structures with adequate fields.
% The following code creates a structure D1 for a drift space
% and a structure QF for a quadrupole.

D1 = atdrift('DR01', 3.0);
QF = atquadrupole('QF', 1.0, 0.2);

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

fodocell = {QF D1 QD D2 QF};

whos fodocell
% LENGTH is useful to find the number of elements in a sequence

L = length(fodocell) 
pause % Press any key to continue;
clc
% Use {:} cell array syntax to print some or all elements
fodocell{1}
pause % FODOCELL{:} will print a long list of all elements. Press any key to continue
clc
fodocell{:}
pause % Press any key to continue;
clc
% Let's build a cell array FODORING that represents a closed ring 
% with 10 periods of FODOCELL the same way we would build 
% any other array in MATLAB from the comman dline

fodoring = [fodocell fodocell fodocell fodocell fodocell...
           fodocell fodocell fodocell fodocell fodocell]; 
        
whos fodoring
pause % Press any key to continue;
clc
% The first element in THERING is 
fodoring{1}

% To inspect or change the value of a specific field we can use MATLAB syntax
% for accessing cells in cell arrays and field in structures
oldK = fodoring{1}.K

fodoring{1}.K = 0.25;

newK = fodoring{1}.K

pause % Press any key to continue;
clc
% Lattice THERING is a variable in MATLAB workspace.
% We can use it in accelerator physics functions and scripts
%
% For example: function FindM44 finds 4-by-4 transverse transfer matrix
M = findm44(fodoring,0)
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
