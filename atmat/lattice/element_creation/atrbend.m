function elem=atrbend(fname,varargin)
%ATRBEND Creates a rectangular bending magnet element with class 'Bend'
%
%  Two calling methods (that can be combined)
%  ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
%
%  INPUTS
%  1. FNAME        - Family name 
%  2. LENGTH       - Length of the arc for an on-energy particle
%                     [m], default to 0
%  3. BENDINGANGLE - Total bending angle [rad], defaults to 0 
%  4. K			   - Focusing strength, defaults to 0
%  5. PASSMETHOD   -Tracking function, defaults to 'BendLinearPass'
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
%  1. ELEM - Structure with the AT element
%
%  EXAMPLES
%  1. atrbend(famname,length,bendingangle,k,passmethod,'fieldname1',value1,...)
%    each pair {'fieldname',value} is added to the element
%
%  NOTES
%  1. Fieldname can be called by calling the passmethod
%     [req opt] = BndMPoleSymplectic4Pass
%                 where req are mandatory field and opt are optional fields
%  2. Model for BndMPoleSymplectic4Pass (Rad) can be selected with extra
%            fields
%
%       FringeBendEntrance/FringeBendExit = 0,1,2,3
%       Version 0 no dipole fringe fields
%       Version 1 legacy version Brown First Order (K. Brown. A First and Second Order 
%                  Matrix Theory for the Design of Beam Transport Systems and Charged 
%                  Particle Spectrometers. Internal report, SLAC-75, 1982)
%       Version 2 SOLEIL close to second order of Brown (J. Bengtsson and M. Meddahi. 
%                 Modeling of Beam Dynamics and Comparison with Measurements for 
%                 the Advanced Light Source. London, UK, 1994.)
%       Version 3 THOMX (Dipole Fringe Field Effects in the ThomX Ring, J. Zhang and 
%                 A. Loulergue, Proceedings of IPAC2013, Shanghai, China)
%
%       FringeQuadEntrance/FringeQuadExit = 0,1,2
%       Version 0 no quadrupole fringe fields
%       Version 1 Lee-Whiting Formula
%       Version 2 Linear quadrupole fringe field using the 5 integrant a la
%                 Elegant          
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

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
