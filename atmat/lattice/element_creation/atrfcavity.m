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
[rsrc,L]=getatarg(rsrc,L,'Length');
[rsrc,V]=getatarg(rsrc,V,'Voltage');
[rsrc,F]=getatarg(rsrc,F,'Frequency');
[rsrc,H]=getatarg(rsrc,H,'HarmNumber');
[rsrc,E]=getatarg(rsrc,E,'Energy');
[rsrc,method]=getatarg(rsrc,method,'PassMethod');
[rsrc,timelag]=getatarg(rsrc,0,'TimeLag');
elem=atbaselem(fname,method,'Class','RFCavity','Length',L,...
    'Voltage',V,'Frequency',F,'HarmNumber',H,'TimeLag',timelag,...
    'Energy',E,rsrc{:});
end
