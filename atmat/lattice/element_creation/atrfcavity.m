function elem=atrfcavity(fname,varargin)
%ATRFCAVITY - creates an rfcavity element with Class 'RFCavity'
%
%  ATRFCAVITY(FAMNAME,LENGTH,VOLTAGE,FREQUENCY,HARMONICNUMBER,ENERGY,PASSMETHOD)
%	
%  INPUTS
%   1. FAMNAME	    Family name
%   2. LENGTH		Length [m]
%   3. VOLTAGE	    Peak voltage [V]
%   4. FREQUENCY	RF frequency [Hz] 
%   5. HARMNUMBER	Harmonic Number
%   6. ENERGY       Energy [eV]
%   7. PASSMETHOD	Tracking function, defaults to 'CavityPass'
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
%    ATRFCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also  ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND, ATSKEWQUAD
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

% Input parser for option
[rsrc,L,V,F,H,E,method] = decodeatargs({0,0,1,1,1.E9,'CavityPass'},varargin);
[L,rsrc]                = getoption(rsrc,'Length',L);
[V,rsrc]                = getoption(rsrc,'Voltage',V);
[F,rsrc]                = getoption(rsrc,'Frequency',F);
[H,rsrc]                = getoption(rsrc,'HarmNumber',H);
[E,rsrc]                = getoption(rsrc,'Energy',E);
[timelag,rsrc]          = getoption(rsrc,'TimeLag',0);
[method,rsrc]           = getoption(rsrc,'PassMethod',method);
[cl,rsrc]               = getoption(rsrc,'Class','RFCavity');

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'Voltage',V,'Frequency',F,'HarmNumber',H,'TimeLag',timelag,...
    'Energy',E,rsrc{:});
end
