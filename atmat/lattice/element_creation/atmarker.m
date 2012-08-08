function Elem=atmarker(fname,method)
%ATMARKER(FAMNAME,PASSMETHOD)
%	creates a marker space element
%
%FAMNAME		family name
%PASSMETHOD     tracking function, defaults to 'IdentityPass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE

if nargin < 2, method='IdentityPass'; end

Elem.FamName=fname;
Elem.PassMethod=method;
