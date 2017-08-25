
function elem = atskewquad(fname,varargin)
%ATSKEWQUAD(FAMNAME,LENGTH,Qs,PASSMETHOD)
%	creates a skew quadrupole element with Class 'Multipole'
%   FAMNAME        family name
%   LENGTH         length
%   Qs              skew quad strength
% atskewquad(Fname, L, Qs, method)
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

 [rsrc,L,Qs,method]=decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[Qs,rsrc]=getoption(rsrc,'Qs',Qs);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[PolynomA,rsrc]=getoption(rsrc,'PolynomA',[0 0]);
[cl,rsrc]=getoption(rsrc,'Class','Multipole');
if ~isempty(Qs), PolynomA(2) = Qs; end

 elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomA',PolynomA,rsrc{:});