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
[PolynomA,rsrc]=getoption(rsrc,'PolynomA',PolynomA);
[PolynomB,rsrc]=getoption(rsrc,'PolynomB',PolynomB);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','ThinMultipole');
elem=atbaselem(fname,method,'Class',cl,'Length',0,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
