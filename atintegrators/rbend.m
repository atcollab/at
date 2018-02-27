function Elem=rbend(fname,L,A,A1,A2,K,method)
%RBEND Creates a rectangular bend in old AT version (Obsolete)
% 
%
%  RBEND('FAMILYNAME',  Length[m], BendingAngle[rad], EntranceAngle[rad],
%       ExitAngle[rad], K, 'METHOD')
%
%  INPUTS
%  1. fname   - family name
%  2. L       - Length of the arc for an on-energy particle [m]
%  3. A		  - Total bending angle [rad]
%  4. A1	  - Entrance angle in [rad] (A/2 - for rectangular bends)
%  5. A2	  - Exit angle in[rad] (A/2 - for rectangular bends)
%  6. K		  - quadrupole K-value for combined funtion bends
%  7. method  - name of the function to use for tracking
% 
%  OUTPUTS
%  1. Elem - returns assigned address in the FAMLIST that is uniquely identifies
%             the family
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
%  See also rbend2, rbend3, atrbend, atsbend

ElemData               = atrbend(fname,L,A,K,method);
ElemData.EntranceAngle = A1;  %for backwards compatibility
ElemData.ExitAngle     = A2;

global FAMLIST
Elem = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{Elem}.FamName = fname;
FAMLIST{Elem}.NumKids = 0;
FAMLIST{Elem}.KidsList= [];
FAMLIST{Elem}.ElemData= ElemData;

