function Elem=rbend3(fname,L,A,A1,A2,K, gap, FF1, FF2, method)
%RBEND3 creates a new family in the FAMLIST - a structure with fields
%
%  INPUTS
%    1. fname   - family name
%	 2. L       - Length of the arc for an on-energy particle [m]
%	 3. A		- Total bending angle [rad]
%	 4. A1	    - Entrance angle in [rad] (A/2 - for rectangular bends)
%	 5. A2		- Exit angle in[rad] (A/2 - for rectangular bends)
%	 6. K		- Quadrupole K-value for combined funtion bends
%    7. FF1     - Entrance fringe field
%    7. FF2     - Exit fringe field
%	 8. gap     - Dipole fringe field
%	 9. method  - Name of the function to use for tracking
% 
%   OUTPUTS
%     1. Elem - Returns assigned address in the FAMLIST that is uniquely identifies
%        the family
%
%   NOTES
%     1. Deprecated function, use atrbend instead
%
%   See also rbend, rbend3, atrbend, atsbend

% 
% Added by Laurent S. Nadolski, SOLEIL, 03/04

ElemData                = atrbend(fname,L,A,K,method);
ElemData.EntranceAngle  = A1;  %for backwards compatibility
ElemData.ExitAngle      = A2;
ElemData.FullGap   		= gap;
ElemData.FringeInt1	    = 0.5*FF1; % same convention as in Tracy II
ElemData.FringeInt2	    = 0.5*FF2; % same convention as in Tracy II

global FAMLIST
Elem = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{Elem}.FamName  = fname;
FAMLIST{Elem}.NumKids  = 0;
FAMLIST{Elem}.KidsList = [];
FAMLIST{Elem}.ElemData = ElemData;
