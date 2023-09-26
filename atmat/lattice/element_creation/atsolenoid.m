function Elem=atsolenoid(fname,varargin)
%ATSOLENOID Creates a new solenoid element with Class 'Solenoid'
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
%  1. ELEM - Structure with the AT element
%
%  NOTES
%  1. Fieldname can be called by calling the passmethod
%     [req opt] = BndMPoleSymplectic4Pass
%                 where req are mandatory field and opt are optional
%                 fields
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atrbend atskewquad,
%          atthinmultipole, atmarker, atcorrector

% Input parser for option
[rsrc,L,K,method]   = decodeatargs({0,0,'SolenoidLinearPass'},varargin);
[L,rsrc]            = getoption(rsrc,'Length',L);
[K,rsrc]            = getoption(rsrc,'KS',K);  % Kept for compatibilty
[K,rsrc]            = getoption(rsrc,'K',K);   % Correct attribute name
[method,rsrc]       = getoption(rsrc,'PassMethod',method);
[cl,rsrc]           = getoption(rsrc,'Class','Solenoid');

% Gradient setting if not specified explicitly
Elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'K',K,rsrc{:});
end
