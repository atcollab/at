function Elem=atsextupole(fname,L,S,method)
%ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD)
%	creates a sextupole element
%
%FAMNAME		family name
%LENGTH			length [m]
%S				strength [m-2]
%PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATMULTIPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 4, method='StrMPoleSymplectic4Pass'; end

Elem.FamName=fname;
Elem.Length=L;
Elem.PolynomA=[0 0 0 0];	 
Elem.PolynomB=[0 0 S 0]; 
Elem.MaxOrder=3;
Elem.NumIntSteps=10;
Elem.PassMethod=method;
