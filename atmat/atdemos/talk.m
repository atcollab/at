%TALK

% keeping notes and comments in lattice files 
help spear2
help spear2rad
help sp3v81loco

% STUDY 1: load a file
spear2
whos

THERING{1:10}
FAMLIST{1:10}

% Read MAD output
clear all
THERING = readmad('spear3_mad_structure.out');



% STUDY 2
% LATTICE MANIPULATION 
% start over by clearing the workspace and loading spear2
clear all
spear2;

% You may access/modify parameters of individual elements
% in several ways
% 1. Text mode is very useful to write your own scripts:
%    for example
THERING{3} % notice curly {} braces for cell access
%	displays the data from the first Q3 quadrupole
% THERING{i} is actually a MATLAB structure so you can 
% access individual fields using "."
THERING{3}.R1 
% displays only the 6x6 entrance tilt matrix R1 for the same quadrupole
% You may even access individual elements in that matrix
THERING{3}.R1(1,1)
% Try
getfield(THERING{3},'PassMethod') % WATCH OUT -  you will get 
% MATLAB warnings if the field does not exist for that element
% or if you misspell it. MATLAB is case sensitive

% SEE ALSO MATLAB help on SETFIELD, GETFIELD, FIELDNAMES   funtions
% useful for structure access / manipulation
fieldnames(THERING{3})
% returns the list of fields in for THERING{3} 
% In fact this list is common for all Q3 



% GUI Lattice manipulation

clear all
spear2 
% Draw the ring - select and edit elements
intlat

% GUI editor of an elemet 
intelem(3)

% GUI (slider) control of a 'K' parameter in
% all QF and QD simultaneously
demoknob


% STUDY3: TRACKING
spear2
% Use RINGPASS for tracking
help ringpass
% Input vectors must be Nx6 matrixes
X0 = [0.01 0 0 0 0 0]'
% Track 1 particle - 1 turn
X = ringpass(THERING,X0);

% Track 1 particle - 10 turns
X = ringpass(THERING,X0,10)

% Track 3 particles - 10 turns
X0 = [X0 X0 X0]
X = ringpass(THERING,X0,10)



% Using results of tracking for accelearator physics calculations in MATLAB
X0 = [0.001 0 0 0 0 0]';
% Track 1 particle - 200 turns
clear X
X = ringpass(THERING,X0,200);
plot(X(1,:),'.')
xlabel('Turn Number');
ylabel('X [mm]');
figure
plot(X(1,:),X(2,:),'.')
title('Phase Space');
figure
xspectrum = abs(fft(X(1,:)));
plot(xspectrum)
title('FFT');


% Plot phase-space near 1/3 integer resonance
clear all; close all
spear2
fittune2([0.345 0.26],'QF','QD');
X0 = [0.001 0 0 0 0 0]'*(1:15);
% Track 15 particles - 200 turns
clear X
X = ringpass(THERING,X0,200);
figure
plot(X(1,:),X(2,:),'.')
title('Phase Space');



% Low-level element pass method
X0 = [0.001 0 0 0 0 0]';
X1 = QuadLinearPass(THERING{5},X0)
X2 = DriftPass(THERING{6},X1)



% Many physics routines call RINGPASS, LINEPASS
% or one of the element-level mex-functions


% Example 1: CLOSED ORBIT CALCULATION
clear all
spear2

plotcod(THERING,0)
% use intlat to set ByError in one of the bends

plotcod(THERING,0)
% Behold horizontally distorted orbit

clear all
clf
spear2
plotcod(THERING,0)
% use intlat to set T1,T2 in one of the quads
% to T1 = [0.001 0 0.001 0 0 0]; T2= -T1;
plotcod(THERING,0)
% Behold orbit distortion in both planes


% Example 2: Make a script to calculate 
% Momentum Compaction Factor

% Lets calculate momentum compaction factor:
% dL = alpha*L*dP/P

% Find off-energy fixed point 
dP = 0.001
fp = findorbit(THERING,dP)
format long

% Build initial condition vector that starts
% on the fixed point
X0 = fp;
X0(5) = dP;
X0(6) = 0;
X0
% Track X0 for 1 turn
X1 = ringpass(THERING,X0)
% Calculate momentum compaction factor
alpha = X1(6)/(dP*234)

% EXAMPLE 3. Beta-functions and tunes 
clear all
spear2
plotbeta
% Use demoknob to tweek quad values
demoknob
plotbeta


