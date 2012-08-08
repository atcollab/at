function Elem=atthinmultipole(fname,PolynomA,PolynomB,method)
%ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD)
%	creates a thin multipole element
%
%	FAMNAME			family name
%	POLYNOMA        skew [dipole quad sext oct];	 
%	POLYNOMB        normal [dipole quad sext oct]; 
%	PASSMETHOD      tracking function. Defaults to 'ThinMPolePass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 4, method='ThinMPolePass'; end
la=length(PolynomA);
lb=length(PolynomB);
lg=max(la,lb);

Elem.FamName=fname;
Elem.PolynomA=[PolynomA zeros(lg-la,1)];	 
Elem.PolynomB=[PolynomB zeros(lg-lb,1)];
Elem.MaxOrder=lg-1;
Elem.PassMethod=method;
