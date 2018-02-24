function elem=atquadrupole(fname,varargin)
%ATQUADRUPOLE - creates a quadrupole element with Class 'Quadrupole'
%
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)	
%
%  INPUTS
%    1. FAMNAME	    family name
%    2. LENGTH      length [m]
%    3. K           strength [m-2]
%    4. PASSMETHOD  tracking function, defaults to 'QuadLinearPass'
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
%    1. ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%       Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND, ATSKEWQUAD,
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR, ATRINGPARAM

[rsrc,L,K,method] = decodeatargs({0,[],'QuadLinearPass'},varargin);
[L,rsrc]          = getoption(rsrc,'Length',L);
[K,rsrc]          = getoption(rsrc,'K',K);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]   = getoption(rsrc,'PolynomB',[0 0]);
[cl,rsrc]         = getoption(rsrc,'Class','Quadrupole');

% Gradient setting if not specified explicitly
if ~isempty(K), PolynomB(2)=K; end

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,'K',PolynomB(2),...
    'PolynomB',PolynomB,rsrc{:});
end
