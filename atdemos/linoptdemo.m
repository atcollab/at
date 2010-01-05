%LINOPTDEMO script illustrates the use of LINOPT
% It reproduces plots of couplig parameters - figures 1,2 for Ref[1]
%
% [1] A.Terebilo, Accelerator Modeling with MATLAB Accelerator Toolbox
%     Proceedings of PAC 2001
%
% See also LINOPT


% The lattice for this type of analysis must NOT contain time-dependent
% elements (RF cavity). Elements should NOT use pass-methods with radiation
% (StrMPoleSymplectic4Pass)
% load spear2 lattice
spear2;

% Find  indexes of elements in QF quadrupole family
QFI = findcells(THERING,'FamName','QF');
% Select some of them to randomly tilt
TILTI = QFI([3:7 10:12]);

% NOTE: How to introduce random coupling and misalignment errors:
% s-rotations(tilts) and transverse displacements (shifts)

% 1.generate random  rotations 
tilterr = 1*pi/180;			% RMS tilt error [degrees]
qftilts = tilterr*randn(1,length(TILTI));

% 2. rotate elements
settilt(TILTI,qftilts);

NE = length(THERING)+1;
LinOptOutput = linopt(THERING,0,1:NE);

% copy LINOPT output 'LinOptOutput' into separate arrays for plotting
GG = cat(1,LinOptOutput.gamma);
spos = cat(1,LinOptOutput.SPos);
CC = reshape(cat(2,LinOptOutput.C),4,[]);

subplot(2,1,1)
plot(spos,CC(1,:),'.-r')
hold on
plot(spos,CC(2,:),'.-k')
plot(spos,CC(3,:),'.-b')
plot(spos,CC(4,:),'.-g')

%Scale axis
SCALE = axis;
SCALE(2) = 0;
SCALE(2) = spos(end);
axis(SCALE);
%Annotate
legend('C_1_1','C_2_1','C_1_2','C_2_2')
title('Elements of coupling matrix C_i_j');
hold off

% Second subplot
subplot(2,1,2)
plot(spos,GG,'.-k')
title('Mixing parameter \gamma')

%Scale 
SCALE = axis;
SCALE(1) = 0;
SCALE(2) = spos(end);
axis(SCALE);

%Annotate
legend('\gamma')
xlabel('s - position [m]')
