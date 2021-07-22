function elem=atquadrupole(fname,varargin)
%ATQUADRUPOLE Creates a quadrupole element with Class 'Quadrupole'
%
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)	
%
%  INPUTS
%  1. FAMNAME    - Family name
%  2. LENGTH     - Length [m]
%  3. K          - Strength [m-2]
%  4. PASSMETHOD - Tracking function, defaults to 'StrMPoleSymplectic4Pass'
%
%  OPTIONS (order does not matter)
%    R1			 -	6 x 6 rotation matrix at the entrance
%	 R2        	 -	6 x 6 rotation matrix at the entrance
%	 T1			 -	6 x 1 translation at entrance 
%	 T2			 -	6 x 1 translation at exit
%	 NumIntSteps -   Number of integration steps
%	 MaxOrder    -   Max Order for multipole (1 up to quadrupole)
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  EXAMPLES
%  1. Fieldname can be called by calling the passmethod
%     [req opt] = StrMPoleSymplectic4Pass
%                 where req are mandatory field and opt are optional fields
%  2. atquadrupole(famname,length,k,passmethod,'fieldname1',value1,...)
%       each pair {'fieldname',value} is added to the element
%
%  3. Quadrupole fringe field can be activated at element entrance or exit
%     with option FringeQuadEntrance/FringeQuadExit=0,1,2
%     Version 0: no fringe field
%     Version 1: Lee-Whiting formula
%     Version 2: Lee-Whiting Elegant-like formula where 5 integral need to
%     be provided
%     
%  See also atdrift, atsextupole, atsbend, atrbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector, atringparam

[rsrc,L,K,method] = decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]          = getoption(rsrc,'Length',L);
[K,rsrc]          = getoption(rsrc,'K',K);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]   = getoption(rsrc,'PolynomB',[0 0]);
[PolynomA,rsrc]   = getoption(rsrc,'PolynomA',[0 0]);
[cl,rsrc]         = getoption(rsrc,'Class','Quadrupole');

% Gradient setting if not specified explicitly
if ~isempty(K), PolynomB(2)=K; end

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,'K',PolynomB(2),...
    'PolynomB',PolynomB,'PolynomA',PolynomA,'DefaultMaxOrder',1,rsrc{:});
end
