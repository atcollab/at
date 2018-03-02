function z = skewquad(fname,L , Qs , method)
%SKEWQUAD Creates a skewquad element (alias to multipole) in old AT version (Obsolete)
%
%skewquad(Fname, L, Qs, method)
%  INPUTS
%  1. FAMNAME - Family name
%  2. LENGTH  - Length in meters [m]
%  3. Qs      - Skew quad strength [m-2]
%  4. method  - Integration passsmethod
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
%  1.  atskewquad(Fname, L, Qs, method)
%
%  NOTES
%  1. Obsolete: consider using atskewquad instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atrbend,
%          atmultipole, atthinmultipole, atmarker, atcorrector atskewquad


z = multipole(fname,L, [0 Qs 0 0], [0 0 0 0], method);