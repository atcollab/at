function elem=atmarker(fname,varargin)
%ATMARKER - creates a marker space element
%
%  ATMARKER(FAMNAME,PASSMETHOD)
%	
%  INPUTS  
%    1. FAMNAME		family name
%    2. PASSMETHOD     tracking function, defaults to 'IdentityPass'
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
%    ATMARKER(FAMNAME,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  NOTES
%      1. Fieldname can be called by calling the passmethod
%         [req opt] = IdentityPass
%                     where req are mandatory field and opt are optional
%                     fields
%
%See also  ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND ATSKEWQUAD,
%          ATTHINMULTIPOLE, ATCORRECTOR

% Input parser for option
[rsrc,method,~] = decodeatargs({'IdentityPass',''},varargin);
[method,rsrc]   = getoption(rsrc,'PassMethod',method);
[cl,rsrc]       = getoption(rsrc,'Class','Marker');

% Build the element
elem=atbaselem(fname,method,'Class',cl,rsrc{:});
end
