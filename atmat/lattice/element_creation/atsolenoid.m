function Elem=atsolenoid(fname,varargin)
%ATSOLENOID - creates a new solenoid element with Class 'Solenoid'
%
%   Elem =solenoid('FAMILYNAME',Length [m],KS,'METHOD')
%	
%  INPUTS
%	1. FamName		  family name
%	2. Length	      length[m]
%	3. KS             solenoid strength KS [rad/m]
%	4. PassMethod     name of the function to use for tracking
%
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
%  NOTES
%      1. Fieldname can be called by calling the passmethod
%         [req opt] = BndMPoleSymplectic4Pass
%                     where req are mandatory field and opt are optional
%                     fields
%
%  See Also ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND ATSKEWQUAD,
%          ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

% Input parser for option
[rsrc,L,KS,method]  = decodeatargs({0,0,'SolenoidLinearPass'},varargin);
[L,rsrc]            = getoption(rsrc,'Length',L);
[KS,rsrc]           = getoption(rsrc,'KS',KS);
[method,rsrc]       = getoption(rsrc,'PassMethod',method);
[cl,rsrc]           = getoption(rsrc,'Class','Solenoid');

% Gradient setting if not specified explicitly
Elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'KS',KS,rsrc{:});
end
