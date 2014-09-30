function elem=atdrift(fname,varargin)
%ATDRIFT(FAMNAME,LENGTH,PASSMETHOD)
%	creates a drift space element with Class 'Drift'
%
%FAMNAME		family name
%LENGTH			length [m]
%PASSMETHOD     tracking function, defaults to 'DriftPass'
%
%ATDRIFT(FAMNAME,LENGTH,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,method]=decodeatargs({0,'DriftPass'},varargin);
[rsrc,L]=getatarg(rsrc,L,'Length');
[rsrc,method]=getatarg(rsrc,method,'PassMethod');
elem=atbaselem(fname,method,'Class','Drift','Length',L,rsrc{:});
end
