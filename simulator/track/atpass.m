%ATPASS is a numerical tracking engine for AT 1.3
% TRJ = ATPASS(LATTICE, RIN, MODE, NTURNS, REFPTS)
%   LATTICE     AT lattice 
%   MODE        0 - reuse lattice
%               1 - new lattice
%   NTURNS      number of turns
%   REFPTS      Indexes of elements where the trajectory is observed
%               May run from 1 to length(LATTICE)+1
%               Ignored if NTURNS > 1
% See also  RINGPASS, LINEPASS