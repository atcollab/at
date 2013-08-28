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
[rsrc,cl]=getatarg(rsrc,'Marker','Class');
elem=atbaselem(fname,method,'Class',cl,'Length',0,rsrc{:});
end
