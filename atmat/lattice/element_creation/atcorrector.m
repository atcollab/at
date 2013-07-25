function Elem=atcorrector(fname,L,kick,method)
%ATCORRECTOR(FAMNAME,LENGTH,KICK,PASSMETHOD)
%	creates a drift space element with class 'Corrector'
%
%FAMNAME		family name
%LENGTH			length [m]
%KICK           [hor. kick vert. kick] [rad]
%PASSMETHOD     tracking function, defaults to 'CorrectorPass'
%
%See also: ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER

if nargin < 4, method='CorrectorPass'; end

Elem.FamName=fname;  % add check for existing identical family names
Elem.Length=L;
Elem.KickAngle=kick;
Elem.PassMethod=method;
Elem.Class='Corrector';
