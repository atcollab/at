function elem=atmultipole(fname,varargin)
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

[rsrc,L,PolynomA,PolynomB,method]=decodeatargs({0,0,0,'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[PolynomA,rsrc]=getoption(rsrc,'PolynomA',PolynomA);
[PolynomB,rsrc]=getoption(rsrc,'PolynomB',PolynomB);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Multipole');
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
