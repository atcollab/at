function elem=atsextupole(fname,varargin)
%ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD)
%	creates a sextupole element with class 'Sextupole'
%
%FAMNAME		family name
%LENGTH			length [m]
%S				strength [m-2]
%PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
%
%ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATQUADRUPOLE, ATMULTIPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,S,method]=decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]=getoption(rsrc,'PolynomB',[0 0 0]);
[cl,rsrc]=getoption(rsrc,'Class','Sextupole');
if ~isempty(S), PolynomB(3)=S; end
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomB',PolynomB,rsrc{:});
end
