function Elem=atquadrupole(fname,L,K,method)
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)
%	creates a quadrupole element
%
%FAMNAME		family name
%LENGTH			length [m]
%K				strength [m-2]
%PASSMETHOD     tracking function, defaults to 'QuadLinearPass'
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 4, method='QuadLinearPass'; end

Elem.FamName=fname;  % add check for existing identical family names
Elem.Length=L;
Elem.PolynomA=[0 0 0];	 
Elem.PolynomB=[0 K 0];
Elem.MaxOrder=2;
Elem.NumIntSteps=10;
Elem.K=K;
Elem.PassMethod=method;
