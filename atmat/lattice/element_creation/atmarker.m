function elem=atmarker(fname,varargin)
%ATMARKER Creates a marker space element
%
%  ATMARKER(FAMNAME,PASSMETHOD)
%	
%  INPUTS  
%  1. FAMNAME	 - Family name
%  2. PASSMETHOD - Tracking function, defaults to 'IdentityPass'
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
%  1. atmarker(famname,passmethod,'fieldname1',value1,...)
%    each pair {'fieldname',value} is added to the element
%
%  NOTES
%  1. Fieldname can be called by calling the passmethod
%     [req opt] = IdentityPass
%                 where req are mandatory field and opt are optional fields
%
%See also  atdrift, atquadrupole, atsextupole, atsbend, atrbend atskewquad,
%          atthinmultipole, atcorrector

% Input parser for option
[rsrc,method,~] = decodeatargs({'IdentityPass',''},varargin);
[method,rsrc]   = getoption(rsrc,'PassMethod',method);
[cl,rsrc]       = getoption(rsrc,'Class','Marker');

% Build the element
elem=atbaselem(fname,method,'Class',cl,rsrc{:});
end
