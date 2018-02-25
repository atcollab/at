function elem=atmonitor(fname,varargin)
%ATMONITOR Creates a Beam Position Monitor element with Class 'Monitor'
%
%  INPUTS
%  1. fname - Family name
% 
%  ATMONITOR(FAMNAME,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also atdrift, atsextupole, atsbend, atrbend
%          atmultipole, atthinmultipole, atmarker, atcorrector

[rsrc,method,~]=decodeatargs({'IdentityPass',''},varargin);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Monitor');
elem=atbaselem(fname,method,'Class',cl,rsrc{:});
end
