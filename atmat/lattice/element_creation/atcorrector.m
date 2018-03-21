function elem=atcorrector(fname,varargin)
%ATCORRECTOR Creates a drift space element with class 'Corrector'
%
%  atcorrector(FAMNAME,LENGTH,KICK,PASSMETHOD)
%	
%  INPUTS
%  1. FAMNAME		family name
%  2. LENGTH		length [m]
%  3. KICK        [hor. kick, vert. kick] [rad]
%  4. PASSMETHOD  tracking function, defaults to 'CorrectorPass'
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  EXAMPLES
%  1. Each pair {'FIELDNAME',VALUE} is added to the element
%
%  NOTES
%  1. Fieldname can be called by calling the passmethod
%     [req opt] = CorrectorPass
%                 where req are mandatory field and opt are optional fields
%
%  See also atquadrupole, atsextupole, atsbend, atrbend
%           atmultipole, atthinmultipole, atmarker

% Input parser for option
[rsrc,L,kick,method] = decodeatargs({0,[0 0],'CorrectorPass'},varargin);
[L,rsrc]             = getoption(rsrc,'Length',L);
[kick,rsrc]          = getoption(rsrc,'KickAngle',kick);
[method,rsrc]        = getoption(rsrc,'PassMethod',method);
[cl,rsrc]            = getoption(rsrc,'Class','Corrector');

% Build the element
elem = atbaselem(fname,method,'Class',cl,'Length',L,'KickAngle',kick,rsrc{:});
end
