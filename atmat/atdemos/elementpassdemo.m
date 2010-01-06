%ELEMENTPASSDEMO self-running tutorial
% 1. Phase-Space tracking variables
% 2. Tracking through individual elements
% 3. Method - Element consistencyclear all
clc
echo on
% The term 'tracking' in Accelerator Physics refers to numerical simulation
% of particle motion in phase-space as it passes through an accelerator 
%   
% MATLAB Accelerator Toolbox uses 6-by-1 column vectors to represent
%  Individual particles in phase space with components [x px y py delta ct]'
%  For example: 
R = [0.01 0 0.01 0 0 0]'

% 6-by-N matrixes are used to represent groups of N particles
% 
RRR = [R 2*R 3*R]

pause % Press any key to continue 
clc
% In Accelerator Toolbox tracking is built upon an open-ended 
% collection of functions that track particles through
% individual accelerator elements
%
% Examle: Load spear2 lattice
spear2
%
% Second element in spear2 lattice is a drift space
SOMEDRIFT = THERING{2}
whos SOMEDRIFT
% SOMEDRIFT is a MATLAB structure
% Now use function DRIFTPASS to track through SOMEDRIFT
pause % Press any key to continue 
clc
driftpass(SOMEDRIFT,R)
%
% DRIFTPASS and other tracking functions in accelerator Toolbox 
% accept matrix input to simultaneously track many particles
% 
driftpass(SOMEDRIFT,RRR)  
% Obviously in a drift space particle momentums don't change
%
% Try this
DriftPass(SOMEDRIFT,[0 0.01 0 0.02 0 0]'),

pause % Press any key to continue 
clc
% Accelerator Toolbox provides an open endeed collection
% of functions that track through elements using various
% field models.
%
% For example with a more interesting element QUADRUPOLE 
% the user can  use different methods
% implemented as different pass-functions:

SOMEQUAD = THERING{5};
% ______________________________________________________
QuadLinearPass(SOMEQUAD,R)
% ______________________________________________________
StrMPoleSymplectic4Pass(SOMEQUAD,R)
% ______________________________________________________
StrMPoleSymplectic4RadPass(SOMEQUAD,R)

pause % Press any key to continue 
clc
% The choice of a proper model depends on
%
% 1. The problem 
%          
% 2. Speed-accuracy trade-off 
%      For example:
%      StrMPoleSymplectic4Pass (4-th order integrator)
%      is slower but more accurate
%      than StrMPoleSymplectic2Pass (2-nd order integrator)     
% 3. Physical considerations 
%      For example:
%      DRIFTPASS assumes a field-free region which is
%      NOT a good model for a quadrupole magnet
% 4. Element-Method consistency
%      Element data gets passed to a pass-function as the first argument
%      Pass-function attempts to use the field with specific name:
%      For example:
%      QUADLINEARPASS needs fields 'Length' and 'K' ...
%      If the element is a drift it does not have 'K' field
%      If in the above examples we tried QUADLINEARPASS(SOMEDRIFT,R)
%      MATLAB would ungracefully stop excecution
%      !!! This feature puts responsibility for consistency between
%      Pass-functions used and elements ON THE USER. Small price to 
%      pay for flexibility !!!
%
pause % Press any key to continue
clc
% Available and extensively tested methods in Accelerator Toolbox 1.0 
%
% AperturePass
% BendLinearPass
% BndMPoleSymplectic4Pass
% BndMPoleSymplectic4RadPass
% DriftPass
% IdentityPass
% QuadLinearPass        
% StrMPoleSymplectic4Pass
% StrMPoleSymplectic4RadPass
% ThinCavityPass
% ThinCorrectorPass
% 
% The names were ment to be long and self-explanatory and end with 'Pass'
%
% Calling syntax is allways for all element pass-functions is the same
%
% These files are originally written in C and converted to MATLAB mex-functions
% They are located (together with source codes and some with help files)
% in ..\simulator\element
                        
pause % Press any key to finish 
clc

echo off
clc