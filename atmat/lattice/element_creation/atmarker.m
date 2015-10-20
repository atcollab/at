function elem=atmarker(fname,varargin)
%ATMARKER(FAMNAME,PASSMETHOD)
%	creates a marker space element
%
%FAMNAME		family name
%PASSMETHOD     tracking function, defaults to 'IdentityPass'
%
%ATMARKER(FAMNAME,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE


[rsrc,method,~]=decodeatargs({'IdentityPass',''},varargin);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Marker');
elem=atbaselem(fname,method,'Class',cl,rsrc{:});
end
