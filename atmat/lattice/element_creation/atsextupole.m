function elem=atsextupole(fname,varargin)
%ATSEXTUPOLE - creates a sextupole element with class 'Sextupole'
%
%  ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD)
%	
%  INPUTS
%	 1. FNAME        	family name 
%    2. LENGTH			length [m]
%    3. S				strength [m-2]
%    4. PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
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
%      1. ELEM - Structure with the AT element
%
%  EXAMPLES
%    ATSEXTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also: ATDRIFT, ATQUADRUPOLE, ATMULTIPOLE, ATSBEND,
%            ATRBEND,ATSKEWQUAD, ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, 
%            ATCORRECTOR

% Input parser for option
[rsrc,L,S,method] = decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]          = getoption(rsrc,'Length',L);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]   = getoption(rsrc,'PolynomB',[0 0 0]);
[cl,rsrc]         = getoption(rsrc,'Class','Sextupole');

% Sextupole setting if not specified explicitly
if ~isempty(S), PolynomB(3)=S; end

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomB',PolynomB,rsrc{:});
end
