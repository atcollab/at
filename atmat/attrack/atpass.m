function varargout=atpass(varargin) %#ok<STOUT>
%ATPASS is a numerical tracking engine for AT 1.3
% ROUT = ATPASS(LATTICE,RIN,MODE,NTURNS,REFPTS)
%
%   LATTICE     AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   MODE        0 - reuse lattice
%               1 - new lattice
%   NTURNS      number of turns
%   REFPTS      Indexes of elements where the trajectory is observed
%               May run from 1 to length(LATTICE)+1
%
% ROUT = ATPASS(LATTICE,RIN,MODE,NTURNS,REFPTS,PREFUNC,POSTFUNC)
%   PREFUNC     function called before each element
%   POSTFUNC    function called after each element
%
% [ROUT,LOST] = ATPASS(LATTICE,RIN,MODE,NTURNS,REFPTS,PREFUNC,POSTFUNC,NHIST)
%   Returns additionally information on lost particles
%
%   NHIST       length of history buffer
%   LOST        1x1 structure with the following fields:
%               turn        1xN vector, turn number where the particle is lost
%               element     1xN vector, element number where the particle is lost
%               coordinates 6xNxNHIST matrix, coordinates at the entrance of the
%                           element where the particle was lost
%
% See also  RINGPASS, LINEPASS
error('at:missingMex','missing MEX file.');
end
