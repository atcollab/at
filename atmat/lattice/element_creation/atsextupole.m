function elem=atsextupole(fname,varargin)
%ATSEXTUPOLE Creates a sextupole element with class 'Sextupole'
%
%  ATSEXTUPOLE(FAMNAME,LENGTH,STRENGTH,PASSMETHOD)
%
%  INPUTS
% 1. FNAME          family name
% 2. LENGTH         length [m]
% 3. STRENGTH       strength [m^{-3}]
% 4. PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
%
%  OPTIONS (order does not matter)
%  R1       6 x 6 rotation matrix at the entrance
%  R2       6 x 6 rotation matrix at the exit
%  T1       6 x 1 translation at entrance
%  T2       6 x 1 translation at exit
%  NumIntSteps    Number of integration steps
%  MaxOrder       Max Order for multipole (2 up to sextupole)
%
%  OUTPUTS
%      1. ELEM - Structure with the AT element
%
%  EXAMPLES
%    ATSEXTUPOLE(FAMNAME,LENGTH,STRENGTH,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also: atdrift, atquadrupole, atmultipole, atsbend,
%            atrbend,atskewquad, atmultipole, atthinmultipole, atmarker,
%            atcorrector, atoctupole

% Input parser for option
[rsrc,L,strength,method] = decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]          = getoption(rsrc,'Length',L);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]   = getoption(rsrc,'PolynomB',[0 0 0]);
[cl,rsrc]         = getoption(rsrc,'Class','Sextupole');

% Sextupole setting if not specified explicitly
if ~isempty(strength), PolynomB(3)=strength; end

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomB','DefaultMaxOrder',2,PolynomB,rsrc{:});
end
