function Elem=atmarker(fname,method)
%ATMARKER(FAMNAME,PASSMETHOD)
%	creates a marker space element with Class 'Marker'
%
%FAMNAME		family name
%PASSMETHOD     tracking function, defaults to 'IdentityPass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE

if nargin < 2, method='IdentityPass'; end

Elem.FamName=fname;
Elem.Length=0;
Elem.PassMethod=method;
Elem.Class='Marker';
