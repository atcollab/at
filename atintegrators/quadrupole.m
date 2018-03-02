function z=quadrupole(fname,L,K,method,varargin)
%QUADRUPOLE Creates a quadrupole element in old AT version (Obsolete)
%quadrupole('familyname',length [m],k,'method')
%  INPUTS
%  1. FAMNAME    - Family name
%  2. LENGTH     - Length [m]
%  3. K          - Strength [m-2]
%  4. PASSMETHOD - Tracking function, defaults to 'QuadLinearPass'
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
%  1. atquadrupole(famname,length,k,passmethod,'fieldname1',value1,...)
%       each pair {'fieldname',value} is added to the element
%
%  NOTES
%  1. Obsolete: consider using atquadrupole instead
%
%  See also atdrift, atsextupole, atsbend, atrbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector, atringparam


ElemData = atquadrupole(fname,L,K,method,varargin{:});

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

