function varargout=atpass(varargin) %#ok<STOUT>
%ATPASS is a numerical tracking engine for AT 1.3
% ROUT = ATPASS(LATTICE, RIN, MODE, NTURNS, REFPTS)
%
%   LATTICE     AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   MODE        0 - reuse lattice
%               1 - new lattice
%   NTURNS      number of turns
%   REFPTS      Indexes of elements where the trajectory is observed
%               May run from 1 to length(LATTICE)+1
%
% [ROUT,LOST] = ATPASS(LATTICE, RIN, MODE, NTURNS, REFPTS)
%   Returns additionally information on lost particles
%
%   LOST        1x1 structure with the following fields:
%               turn        1xN vector, turn number where the particle is lost
%               element     1xN vector, element number where the particle is lost
%               coordinate  6xN matrix, coordinates when the particle was
%                           marked as lost
% See also  RINGPASS, LINEPASS
error('at:missingMex','missing MEX file.');
end
