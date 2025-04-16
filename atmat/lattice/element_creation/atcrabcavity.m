function elem=atcrabcavity(fname,varargin)
%ATCRABCAVITY Creates an crab cavity element with Class 'CrabCavity'
%
%  ATCRABCAVITY(FAMNAME,LENGTH,VOLTAGES,FREQUENCY,HARMONICNUMBER)
%	
%  INPUTS
%   1. FAMNAME	    Family name
%   2. LENGTH		Length [m]
%   3. VOLTAGES	    Array [Horizontal, Vertical] Peak voltages [V]
%   4. FREQUENCY	RF frequency [Hz] 
%   5. HARMNUMBER	Harmonic Number
%
%  OUTPUTS
%      1. ELEM - Structure with the AT element
%
%  EXAMPLES
%    ATCRABCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also  atdrift, atsextupole, atsbend, atrbend, atskewquad, atrfcavity
%          atmultipole, atthinmultipole, atmarker, atcorrector

% Input parser for option
[rsrc,L,V,F,H,method] = decodeatargs({0,0,1,1,'CrabCavityPass'},varargin);
[L,rsrc]                = getoption(rsrc,'Length',L);
[V,rsrc]                = getoption(rsrc,'Voltage',V);
[F,rsrc]                = getoption(rsrc,'Frequency',F);
[H,rsrc]                = getoption(rsrc,'HarmNumber',H);
[timelag,rsrc]          = getoption(rsrc,'TimeLag',0);
[method,rsrc]           = getoption(rsrc,'PassMethod',method);
[cl,rsrc]               = getoption(rsrc,'Class','CrabCavity');

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'Voltages',V,'Frequency',F,'HarmNumber',H,'TimeLag',timelag,...
    rsrc{:});
end
