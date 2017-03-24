function elem=atrbend(fname,varargin)
% ATRBEND creates a rectangular bending magnet element with class 'Bend'
%
%  Two calling methods (that can be combined)
%  ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
%  INPUTS
%	 1. fname        	family name 
%	 2. LENGTH         	length of the arc for an on-energy particle
%                     	[m], default to 0
%	 3. BENDINGANGLE	total bending angle [rad], defaults to 0 
%	 4. K				focusing strength, defaults to 0
%	 5. PASSMETHOD      tracking function, defaults to 'BendLinearPass'
%
%  ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%  Each pair {'FIELDNAME',VALUE} is added to the element
%
%  OUTPUTS
%      1. elem - Structure with the AT element
%
%  NOTES
%      1. Fieldname can be called by calling the passmethod
%         [req opt] = BndMPoleSymplectic4Pass
%                     where req are mandatory field and opt are optional
%                     fields
%
%  See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

% Input parser for option
[rsrc,L,A,K,method]  = decodeatargs({0,0,[],'BendLinearPass'},varargin);
[L,rsrc]             = getoption(rsrc,'Length',L);
[A,rsrc]             = getoption(rsrc,'BendingAngle',A);
[K,rsrc]             = getoption(rsrc,'K',K);
[method,rsrc]        = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]      = getoption(rsrc,'PolynomB',[0 0]);
[cl,rsrc]            = getoption(rsrc,'Class','Bend');
[EntranceAngle,rsrc] = getoption(rsrc,'EntranceAngle',0.5*A);
[ExitAngle,rsrc]     = getoption(rsrc,'ExitAngle',0.5*A);

% Gradient setting if not specified explicitly
if ~isempty(K), PolynomB(2) = K; end

% Build the element
elem = atbaselem(fname,method,'Class',cl,'Length',L,...
    'BendingAngle',A,'EntranceAngle',EntranceAngle,'ExitAngle',ExitAngle,...
    'K',PolynomB(2),'PolynomB',PolynomB,rsrc{:});
end
