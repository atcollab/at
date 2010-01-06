function varargout =  fitchrom2(newchrom, sextfam1, sextfam2, varargin);
%FITCHROM2 fits chromaticity  of THERING using 2 sextupole families
% FITCHROM2(NEWCHROM,SEXTUPOLEFAMILY1,SEXTUPOLEFAMILY2)

% Must declare THERING as global in order for the function to modify sextupole values 
global THERING

%make a column vector
newchrom = newchrom(:);
deltaS = 1e-5; % step size in Sextupole strngth
deltaP = 1e-8;

% find indexes of the 2 quadrupole families use for fitting
S1I = findcells(THERING,'FamName',sextfam1);
S2I = findcells(THERING,'FamName',sextfam2);
InitialS1 = getcellstruct(THERING,'PolynomB',S1I,3);
InitialS2 = getcellstruct(THERING,'PolynomB',S2I,3);

% Compute initial tunes and chromaticities before fitting 

[ LD, InitialTunes] = linopt(THERING,0);
[ LDdP, ITdP] =linopt(THERING,deltaP);

InitialChrom = (ITdP-InitialTunes)/deltaP;

TempTunes = InitialTunes;
TempChrom = InitialChrom;
TempS1 = InitialS1; 
TempS2 = InitialS2;

for i=1:5
		
	% Take Derivative
	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1+deltaS,3);
	[LD , Tunes_dS1 ] = linopt(THERING,0);
	[LD , Tunes_dS1dP ] = linopt(THERING,deltaP);

	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1,3);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2+deltaS,3);
	[LD , Tunes_dS2 ] = linopt(THERING,0);
	[LD , Tunes_dS2dP ] = linopt(THERING,deltaP);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2,3);

	%Construct the Jacobian
	Chrom_dS1 = (Tunes_dS1dP-Tunes_dS1)/deltaP;
	Chrom_dS2 = (Tunes_dS2dP-Tunes_dS2)/deltaP;

	J = ([Chrom_dS1(:) Chrom_dS2(:)] - [TempChrom(:) TempChrom(:)])/deltaS;
	Jinv = inv(J);

	dchrom = (newchrom(:) - TempChrom(:));
	dS = Jinv*dchrom;

	TempS1 = TempS1+dS(1);
	TempS2 = TempS2+dS(2);

	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1,3);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2,3);

	[ LD, TempTunes] = linopt(THERING,0);
	[ LD, TempTunesdP] = linopt(THERING,deltaP);
	TempChrom = (TempTunesdP-TempTunes)/deltaP;
    %disp(TempChrom);

end