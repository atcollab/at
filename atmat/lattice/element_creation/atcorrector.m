function elem=atcorrector(fname,varargin)
%ATCORRECTOR(FAMNAME,LENGTH,KICK,PASSMETHOD)
%	creates a drift space element with class 'Corrector'
%
%FAMNAME		family name
%LENGTH			length [m]
%KICK           [hor. kick vert. kick] [rad]
%PASSMETHOD     tracking function, defaults to 'CorrectorPass'
%
%ATCORRECTOR(FAMNAME,LENGTH,KICK,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER

[rsrc,L,kick,method]=decodeatargs({0,[0 0],'CorrectorPass'},varargin);
elem=atbaselem(fname,method,'Class','Corrector','Length',L,...
    'KickAngle',kick,rsrc{:});
end
