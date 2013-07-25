function Elem=atmonitor(fname)
%ATMONITOR(FAMNAME)
%	creates a Beam Position Monitor element with Class 'Monitor'
%
%FAMNAME		family name
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR



Elem.FamName=fname;  % add check for existing identical family names
Elem.Length=0;
Elem.K=K;
Elem.PassMethod='IdentityPass';
Elem.Class='Monitor';