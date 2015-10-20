function elem=atrfcavity(fname,varargin)
%ATRFCAVITY(FAMNAME,LENGTH,VOLTAGE,FREQUENCY,HARMONICNUMBER,ENERGY,PASSMETHOD)
%	creates an rfcavity element with Class 'RFCavity'
%
%FAMNAME	family name
%LENGTH		length [m]
%VOLTAGE	peak voltage [V]
%FREQUENCY	RF frequency [Hz] 
%HARMNUMBER	Harmonic Number
%ENERGY     Energy [eV]
%PASSMETHOD	tracking function, defaults to 'CavityPass'
%
%ATRFCAVITY(FAMNAME,...,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,V,F,H,E,method]=decodeatargs({0,0,1,1,1.E9,'CavityPass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[V,rsrc]=getoption(rsrc,'Voltage',V);
[F,rsrc]=getoption(rsrc,'Frequency',F);
[H,rsrc]=getoption(rsrc,'HarmNumber',H);
[E,rsrc]=getoption(rsrc,'Energy',E);
[timelag,rsrc]=getoption(rsrc,'TimeLag',0);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','RFCavity');
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'Voltage',V,'Frequency',F,'HarmNumber',H,'TimeLag',timelag,...
    'Energy',E,rsrc{:});
end
