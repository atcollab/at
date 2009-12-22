function [dK, Jinv] =  fittune2(newtunes, quadfam1, quadfam2, varargin);
%FITTUNE2 fits linear tunes of THERING using 2 quadrupole families
% FITTUNE2(NEWTUNES,QUADFAMILY1,QUADFAMILY2)
% [dK, Jinv] = FITTUNE2(NEWTUNES,QUADFAMILY1,QUADFAMILY2)

% Must declare THERING as global in order for the function to modify quadrupole values 
global THERING
if nargin > 3 % use externally supplied step size for quadrupole K-values 
    delta = varargin{1}
else
    delta = 1e-6; % default step size for quadrupole K-values 
end
% find indexes of the 2 quadrupole families use for fitting
Q1I = findcells(THERING,'FamName',quadfam1);
Q2I = findcells(THERING,'FamName',quadfam2);

InitialK1 = getcellstruct(THERING,'K',Q1I);
InitialK2 = getcellstruct(THERING,'K',Q2I);


% Compute initial tunes before fitting 
[ LD, InitialTunes] = linopt(THERING,0);

TempTunes = InitialTunes;
TempK1 = InitialK1;
TempK2 = InitialK2;


% Take Derivative
THERING = setcellstruct(THERING,'K',Q1I,TempK1+delta);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempK1+delta,2);
[LD , Tunes_dK1 ] = linopt(THERING,0);
THERING = setcellstruct(THERING,'K',Q1I,TempK1);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempK1,2);
THERING = setcellstruct(THERING,'K',Q2I,TempK2+delta);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempK2+delta,2);
[LD , Tunes_dK2 ] = linopt(THERING,0);
THERING = setcellstruct(THERING,'K',Q2I,TempK2);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempK2,2);

%Construct the Jacobian
J = ([Tunes_dK1(:) Tunes_dK2(:)] - [TempTunes(:) TempTunes(:)])/delta;
Jinv = inv(J);

dnu = (newtunes(:) - TempTunes(:));
dK = Jinv*dnu;

TempK1 = TempK1+dK(1);
TempK2 = TempK2+dK(2);


THERING = setcellstruct(THERING,'K',Q1I,TempK1);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempK1,2);
THERING = setcellstruct(THERING,'K',Q2I,TempK2);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempK2,2);

[LD,TempTunes] = linopt(THERING,0);
S = sprintf('New tunes = %f %f\n', TempTunes(1), TempTunes(2));
disp(S)

