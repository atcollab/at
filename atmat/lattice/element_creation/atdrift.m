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
[L,rsrc]=getoption(rsrc,'Length',L);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Drift');
elem=atbaselem(fname,method,'Class',cl,'Length',L,rsrc{:});
end
