function elem=atoctupole(fname,varargin)
%ATOCTUPOLE Creates an octupole element with class 'Octupole'
%
%  ATOCTUPOLE(FAMNAME,LENGTH,S,PASSMETHOD)
%	
%  INPUTS
%	 1. FNAME        	family name 
%    2. LENGTH			length [m]
%    3. O				strength [m-3]
%    4. PASSMETHOD     tracking function, defaults to 'StrMPoleSymplectic4Pass'
%
%  OPTIONS (order does not matter)
%    R1				6 x 6 rotation matrix at the entrance
%	 R2        		6 x 6 rotation matrix at the entrance
%	 T1				6 x 1 translation at entrance 
%	 T2				6 x 1 translation at exit
%	 NumIntSteps    Number of integration steps
%	 MaxOrder       Max Order for multipole (1 up to quadrupole)
%
%  OUTPUTS
%      1. ELEM - Structure with the AT element
%
%  EXAMPLES
%    ATOCTUPOLE(FAMNAME,LENGTH,O,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%    Each pair {'FIELDNAME',VALUE} is added to the element
%
%  See also: atdrift, atquadrupole, atsecxtupole, atmultipole, atsbend,
%            atrbend,atskewquad, atmultipole, atthinmultipole, atmarker, 
%            atcorrector

% Input parser for option
[rsrc,L,O,method] = decodeatargs({0,[],'StrMPoleSymplectic4Pass'},varargin);
[L,rsrc]          = getoption(rsrc,'Length',L);
[method,rsrc]     = getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]   = getoption(rsrc,'PolynomB',[0 0 0 0]);
[cl,rsrc]         = getoption(rsrc,'Class','Octupole');

% Octupole setting if not specified explicitly
if ~isempty(O), PolynomB(4)=O; end

% Build the element
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'PolynomB','DefaultMaxOrder',3,PolynomB,rsrc{:});
end
