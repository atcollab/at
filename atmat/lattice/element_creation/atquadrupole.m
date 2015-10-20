function elem=atquadrupole(fname,varargin)
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)
%	creates a quadrupole element with Class 'Quadrupole'
%
%FAMNAME	family name
%LENGTH     length [m]
%K          strength [m-2]
%PASSMETHOD	tracking function, defaults to 'QuadLinearPass'
%
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,K,method]=decodeatargs({0,[],'QuadLinearPass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[K,rsrc]=getoption(rsrc,'K',K);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]=getoption(rsrc,'PolynomB',[0 0]);
[cl,rsrc]=getoption(rsrc,'Class','Quadrupole');
if ~isempty(K), PolynomB(2)=K; end
elem=atbaselem(fname,method,'Class',cl,'Length',L,'K',PolynomB(2),...
    'PolynomB',PolynomB,rsrc{:});
end
