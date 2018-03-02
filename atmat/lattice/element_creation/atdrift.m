function elem=atdrift(fname,varargin)
%ATDRIFT Creates a drift space element with Class 'Drift'
%ATDRIFT(FAMNAME,LENGTH,PASSMETHOD)
%
%  INPUTS
%  1. FAMNAME	   - Family name
%  2. LENGTH	   - Length [m]
%  3. PASSMETHOD - Tracking function, defaults to 'DriftPass'
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
%  1. atdrift(famname,length,passmethod,'fieldname1',value1,...)
%    each pair {'fieldname',value} is added to the element
%

%See also  atquadrupole, atsextupole, atsbend, atrbend
%          atmultipole, atthinmultipole, atmarker, atcorrector

[rsrc,L,method]=decodeatargs({0,'DriftPass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Drift');
elem=atbaselem(fname,method,'Class',cl,'Length',L,rsrc{:});
end
