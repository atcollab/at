function Elem=rbend3(fname,L,A,A1,A2,K, gap, FF1, FF2, method)
%RBEND3 Creates a rectangular bend with different fringe fields at entrance 
%and exit in old AT versions (Obsolete)
%
%  INPUTS
%  1. fname   - family name
%  2. L       - Length of the arc for an on-energy particle [m]
%  3. A		- Total bending angle [rad]
%  4. A1	    - Entrance angle in [rad] (A/2 - for rectangular bends)
%  5. A2		- Exit angle in[rad] (A/2 - for rectangular bends)
%  6. K		- Quadrupole K-value for combined funtion bends
%  7. FF1     - Entrance fringe field
%  8. FF2     - Exit fringe field
%  9. gap     - Dipole fringe field
%  10. method  - Name of the function to use for tracking
% 
%  OUTPUTS
%  1. Elem - Returns assigned address in the FAMLIST that is uniquely identifies
%        the family
%
%  NOTES
%  1. Deprecated function, use atrbend instead
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
%  See also rbend, rbend2, atrbend, atsbend

% 
% Added by Laurent S. Nadolski, SOLEIL, 03/04

ElemData                = atrbend(fname,L,A,K,method);
ElemData.EntranceAngle  = A1;  %for backwards compatibility
ElemData.ExitAngle      = A2;
ElemData.FullGap   		= gap;
ElemData.FringeInt1	    = 0.5*FF1; % same convention as in Tracy II
ElemData.FringeInt2	    = 0.5*FF2; % same convention as in Tracy II

ElemData.FringeBendEntrance	= 2;
ElemData.FringeBendExit 	= 2;
ElemData.FringeQuadEntrance = 0;
ElemData.FringeQuadExit     = 0; 

global FAMLIST
Elem = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{Elem}.FamName  = fname;
FAMLIST{Elem}.NumKids  = 0;
FAMLIST{Elem}.KidsList = [];
FAMLIST{Elem}.ElemData = ElemData;
