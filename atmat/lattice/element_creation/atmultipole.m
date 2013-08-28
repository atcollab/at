function elem=atmultipole(fname,L,varargin)
%ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD)
%	creates a multipole element
%
%	FAMNAME			family name
%	LENGTH			length[m]
%	POLYNOMA        skew [dipole quad sext oct];	 
%	POLYNOMB        normal [dipole quad sext oct]; 
%	PASSMETHOD      tracking function. Defaults to 'StrMPoleSymplectic4Pass'
%
%ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,PolynomA,PolynomB,method]=decodeatargs({0,0,'StrMPoleSymplectic4Pass'},varargin);
[rsrc,PolynomA]=getatarg(rsrc,PolynomA,'PolynomA');
[rsrc,PolynomB]=getatarg(rsrc,PolynomB,'PolynomB');
elem=atbaselem(fname,method,'Class','Multipole','Length',L,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
