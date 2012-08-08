function Elem=atmultipole(fname,L,PolynomA,PolynomB,method)
%ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD)
%	creates a multipole element
%
%	FAMNAME			family name
%	LENGTH			length[m]
%	POLYNOMA        skew [dipole quad sext oct];	 
%	POLYNOMB        normal [dipole quad sext oct]; 
%	PASSMETHOD      tracking function. Defaults to 'StrMPoleSymplectic4Pass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 5, method='StrMPoleSymplectic4Pass'; end
la=length(PolynomA);
lb=length(PolynomB);
lg=max(la,lb);

Elem.FamName=fname;
Elem.Length=L;
Elem.PolynomA=[PolynomA zeros(lg-la,1)];	 
Elem.PolynomB=[PolynomB zeros(lg-lb,1)];
Elem.MaxOrder=lg-1;
Elem.NumIntSteps=10;
Elem.PassMethod=method;
