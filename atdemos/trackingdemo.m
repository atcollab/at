%TRACKINGDEMO self-running tutorial
% 1. Phase-Space tracking variables
% 2. Tracking through individual elements
% 3. Method - Element consistency
clear all
clc
echo on
% The term 'tracking' in Accelerator Physics refers to numerical simulation
% of particle motion in phase-space as it passes through an accelerator 
%   
% MATLAB Accelerator Toolbox uses 6-by-1 column vectors to represent
%  Individual particles in phase space with components [x px y py delta ct]'
%  For example: 
R = [0.01 0 0.01 0 0 0]'

% and 6-by-N matrixes to simultaneously track groups of N particles
%
RRR = [R 2*R 3*R]

pause % Press any key to continue 
clc
% In Accelerator Toolbox tracking is built upon 
% a collection of functions that track particles through
% individual accelerator elements
%
% Example: Load the spear2 lattice
spear2
%
% Second element in spear2 lattice is a drift space
SOMEDRIFT = THERING{2}
whos SOMEDRIFT
% D is a MATLAB structure
% Now use function DRIFTPASS to track through SOMEDRIFT
pause % Press any key to continue 
clc
DriftPass(SOMEDRIFT,R)
% or simultaneously for 3 particles  
DriftPass(SOMEDRIFT,R)  
% Obviously in a drift space particle momentums don't change
%
% Try this
DriftPass(SOMEDRIFT,[0 0.01 0 0.02 0 0]'),

pause % Press any key to continue 
clc
% Accelerator Toolbox provides an open ended collection
% of functions that track through elements using various
% field models.
%
% For example with a more interesting element QUADRUPOLE 
% the user can  use different models
% implemented in as different pass-methds:

SOMEQUAD = THERING{5};
% ______________________________________________________
QuadLinearPass(SOMEQUAD,R)
% ______________________________________________________
StrMPoleSymplectic4Pass(SOMEQUAD,R)
% ______________________________________________________
StrMPoleSymplectic4RadPass(SOMEQUAD,R)
% ______________________________________________________
% even
DriftPass(SOMEQUAD,R)

pause % Press any key to continue 
clc
% The choice of a proper model depends on
%
% 1. The problem 
%          
% 2. Speed-Accuracy trade-off 
%      StrMPoleSymplectic4Pass is slower but more accurate
%      than StrMPoleSymplectic2Pass.      
% 3. Physical considerations 
%      DriftPass assumes a field-free region which is
%      NOT a good model for a quadrupole magnet
% 4. Element-Method consistency
%      Element data gets passed to a pass-function as the first argument
%      Pass-function attempts to use the field with specific name:
%      For example QUADLINEARPASS needs 'Length', 'K', 
%      If the element is a drift it does not have 'K'
%      If in the above examples we tried QUADLINEARPASS(SOMEDRIFT,R)
%      MATLAB would ungracefully stop execution
%      !!! This feature puts responsibility for consistency between
%      Pass-functions used and elements ON THE USER. Small price to 
%      pay for flexibility !!!



pause % Press any key to finish 
echo off
clc
