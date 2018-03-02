function elem=atrfcavity(fname,varargin)
%ATRFCAVITY Creates an rfcavity element with Class 'RFCavity'
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
%  OUTPUTS
%      1. ELEM - Structure with the AT element
%
%  EXAMPLES
%    ATRFCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  NOTES
%      1. Fieldname can be called by calling the passmethod
%         [req opt] = atrfcavity
%                     where req are mandatory field and opt are optional
%                     fields
%%See also  atdrift, atsextupole, atsbend, atrbend, atskewquad
%          atmultipole, atthinmultipole, atmarker, atcorrector

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
