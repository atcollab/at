function z=sbend(fname,L,A,A1,A2,K,method)
%SBEND Creates a sector bend element in old AT versions (Obsolete)
%BEND('FAMILYNAME',  Length[m], BendingAngle[rad], EntranceAngle[rad],
%	ExitAngle[rad], K, 'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName        	family name
%		Length         	length of the arc for an on-energy particle [m]
%		BendingAngle		total bending angle [rad]
%		EntranceAngle		[rad] (0 - for sector bends)
%		ExitAngle			[rad] (0 - for sector bends)
%		ByError				error in the dipole field relative to the design value 
%		K						quadrupole K-value for combined funtion bends
%		PassMethod        name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
%
%  NOTES
%  1. Obsolete: use atsbend instead
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

ElemData=atsbend(fname,L,A,K,method);
ElemData.EntranceAngle=A1;  %for backwards compatibility
ElemData.ExitAngle=A2;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

