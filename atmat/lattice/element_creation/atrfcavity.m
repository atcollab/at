function elem=atrfcavity(fname,L,V,F,H,E,varargin)
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

[rsrc,method,~]=decodeatargs({'CavityPass',''},varargin);
[rsrc,timelag]=getatarg(rsrc,0,'TimeLag');
elem=atbaselem(fname,method,'Class','RFCavity','Length',L,...
    'Voltage',V,'Frequency',F,'HarmNumber',H,'TimeLag',timelag,...
    'Energy',E,rsrc{:});
end
