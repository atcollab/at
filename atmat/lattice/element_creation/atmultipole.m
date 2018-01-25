function elem=atmultipole(fname,varargin)
%ATMULTIPOLE - creates a multipole element
%
%  ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD)
%	
%  INPUTS
%	 1. FNAME        	family name 
%	 2. LENGTH			length[m]
%	 3. POLYNOMA        skew [dipole quad sext oct];	 
%	 4. POLYNOMB        normal [dipole quad sext oct]; 
%	 5. PASSMETHOD      tracking function. Defaults to 'StrMPoleSymplectic4Pass'
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
%    ATMULTIPOLE(FAMNAME,LENGTH,POLYNOMA,POLYNOMB,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND, ATSKEWQUAD,
%          ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

% Input parser for option
[rsrc,L,PolynomA,PolynomB,method] = decodeatargs({0,0,0,'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]                          = getoption(rsrc,'Length',L);
[PolynomA,rsrc]                   = getoption(rsrc,'PolynomA',PolynomA);
[PolynomB,rsrc]                   = getoption(rsrc,'PolynomB',PolynomB);
[method,rsrc]                     = getoption(rsrc,'PassMethod',method);
[cl,rsrc]                         = getoption(rsrc,'Class','Multipole');

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
