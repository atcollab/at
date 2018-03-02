function elem = atskewquad(fname,varargin)
%ATSKEWQUAD Creates a skew quadrupole element with Class 'Multipole'
%atskewquad(famname,length,qs,passmethod)
%
%  INPUTS
%  1. FAMNAME - Family name
%  2. LENGTH  - Length [m]
%  3. Qs      - Skew quad strength [m-2]
%
%  OPTIONS (order does not matter)
%    R1				6 x 6 rotation matrix at the entrance
%	 R2        		6 x 6 rotation matrix at the entrance
%	 T1				6 x 1 translation at entrance 
%	 T2				6 x 1 translation at exit
%	 NumIntSteps    Number of integration steps
%	 MaxOrder       Max Order for multipole (1 up to quadrupole)
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  EXAMPLES 
%  1.  atskewquad(Fname, L, Qs, method)
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atrbend,
%          atmultipole, atthinmultipole, atmarker, atcorrector

% Input parser for option
[rsrc,L,Qs,method] = decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]           = getoption(rsrc,'Length',L);
[Qs,rsrc]          = getoption(rsrc,'Qs',Qs);
[method,rsrc]      = getoption(rsrc,'PassMethod',method);
[PolynomA,rsrc]    = getoption(rsrc,'PolynomA',[0 0]);
[cl,rsrc]          = getoption(rsrc,'Class','Multipole');

% Skew Gradient setting if not specified explicitly
if ~isempty(Qs), PolynomA(2) = Qs; end

% Build the element
 elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomA',PolynomA,rsrc{:});
end