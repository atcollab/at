function elem=atthinmultipole(fname,varargin)
%ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD)
%	creates a thin multipole element
%
%	FAMNAME			family name
%	POLYNOMA        skew [dipole quad sext oct];	 
%	POLYNOMB        normal [dipole quad sext oct]; 
%	PASSMETHOD      tracking function. Defaults to 'ThinMPolePass'
%
%ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,PolynomA,PolynomB,method]=decodeatargs({0,0,'ThinMPolePass'},varargin);
[rsrc,PolynomA]=getatarg(rsrc,PolynomA,'PolynomA');
[rsrc,PolynomB]=getatarg(rsrc,PolynomB,'PolynomB');
[rsrc,method]=getatarg(rsrc,method,'PassMethod');
elem=atbaselem(fname,method,'Class','ThinMultipole','Length',0,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
