function elem=atthinmultipole(fname,varargin)
% ATTHINMULTIPOLE  - creates a thin multipole element
%
% ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD)
%	
%  INPUTS
%	 1. FNAME        	family name 
%	 2. POLYNOMA        skew [dipole quad sext oct];	 
%	 3. POLYNOMB        normal [dipole quad sext oct]; 
%	 4. PASSMETHOD      tracking function. Defaults to 'StrMPoleSymplectic4Pass'
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
%    ATTHINMULTIPOLE(FAMNAME,POLYNOMA,POLYNOMB,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  NOTES
%      1. Fieldname can be called by calling the passmethod
%         [req opt] = BndMPoleSymplectic4Pass
%                     where req are mandatory field and opt are optional
%                     fields
%
%See also  ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND ATSKEWQUAD,
%          ATMULTIPOLE, ATMARKER, ATCORRECTOR

% Input parser for option
[rsrc,PolynomA,PolynomB,method] = decodeatargs({0,0,'ThinMPolePass'},varargin);
[PolynomA,rsrc]                 = getoption(rsrc,'PolynomA',PolynomA);
[PolynomB,rsrc]                 = getoption(rsrc,'PolynomB',PolynomB);
[method,rsrc]                   = getoption(rsrc,'PassMethod',method);
[cl,rsrc]                       = getoption(rsrc,'Class','ThinMultipole');

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',0,...
    'PolynomA',PolynomA,'PolynomB',PolynomB,rsrc{:});
end
