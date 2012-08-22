function Elem=atrfcavity(fname,L,V,F,H,method)
%ATRFCAVITY(FAMNAME,LENGTH,VOLTAGE,FREQUENCY,HARMONICNUMBER,PASSMETHOD)
%	creates an rfcavity element
%
%		FamName			family name
%		Length			length[m]
%		Voltage			peak voltage (V)
%		Frequency		RF frequency [Hz] 
% 		HarmNumber		Harmonic Number
%		PassMethod		name of the function on disk to use for tracking
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 6, method='CavityPass'; end

Elem.FamName = fname;  % add check for identical family names
Elem.Length = L;
Elem.Voltage = V;
Elem.Frequency = F;
Elem.HarmNumber = H; %is this required?  not a required arg in pass method
Elem.PhaseLag = 0;
Elem.PassMethod=method;
Elem.Class='RFCavity';