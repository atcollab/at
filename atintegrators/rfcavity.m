function z = rfcavity(fname,L,V,F,H,method)
%RFCAVITY('FAMILYNAME',Length [m],Voltage[V], Frequency[Hz], Harmonic Number,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length			length[m]
%		Voltage			peak voltage (V)
%		Frequency		RF frequency [Hz] 
% 		HarmNumber		Harmonic Number
%		PassMethod		name of the function on disk to use for tracking
% returns assigned address in the FAMLIST that uniquely identifies
% the family

ElemData = atrfcavity(fname,L,V,F,H,method);

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;


