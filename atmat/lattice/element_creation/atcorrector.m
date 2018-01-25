function elem=atcorrector(fname,varargin)
%ATCORRECTOR - creates a drift space element with class 'Corrector'
%
%  ATCORRECTOR(FAMNAME,LENGTH,KICK,PASSMETHOD)
%	
%  INPUTS
%    1. FAMNAME		family name
%    2. LENGTH			length [m]
%    3. KICK           [hor. kick vert. kick] [rad]
%    4. PASSMETHOD     tracking function, defaults to 'CorrectorPass'
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
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%           ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER

% Input parser for option
[rsrc,L,kick,method] = decodeatargs({0,[0 0],'CorrectorPass'},varargin);
[L,rsrc]             = getoption(rsrc,'Length',L);
[kick,rsrc]          = getoption(rsrc,'KickAngle',kick);
[method,rsrc]        = getoption(rsrc,'PassMethod',method);
[cl,rsrc]            = getoption(rsrc,'Class','Corrector');

% Build the element
elem = atbaselem(fname,method,'Class',cl,'Length',L,'KickAngle',kick,rsrc{:});
end
