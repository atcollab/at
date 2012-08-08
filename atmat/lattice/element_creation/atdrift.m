function Elem=atdrift(fname,L,method)
%ATDRIFT(FAMNAME,LENGTH,PASSMETHOD)
%	creates a drift space element
%
%FAMNAME		family name
%LENGTH			length [m]
%PASSMETHOD     tracking function, defaults to 'DriftPass'
%
%See also: ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 3, method='DriftPass'; end

Elem.FamName=fname;  % add check for existing identical family names
Elem.Length=L;
Elem.PassMethod=method;
