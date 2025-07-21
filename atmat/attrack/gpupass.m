function varargout=gpupass(varargin) %#ok<STOUT>
%GPUPASS Low-level tracking engine on GPU
%
% ROUT = GPUPASS(LATTICE,RIN,MODE,NTURNS,REFPTS,TURN,KEEPCOUNTER,GPUPOOL,INTEGRATOR)
%   LATTICE     AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   MODE        0 - reuse lattice
%               1 - new lattice
%   NTURNS      number of turns
%   REFPTS      Indexes of elements where the trajectory is observed
%               May run from 1 to length(LATTICE)+1
%   KEEPCOUNTER 0 - Start counting turn from 0
%               1 - Keep last turn counter of the previous gpupass call
%   GPUPOOL     GPU to use (see gpuinfo)
%   INTEGRATOR  Type of integrator to use
%               1: Euler 1st order, 1 drift/1 kick per step
%               2: Verlet 2nd order, 1 drift/2 kicks per step
%               3: Ruth 3rd order, 3 drifts/3 kicks per step
%               4: Forest/Ruth 4th order, 4 drifts/3 kicks per step (Default)
%               5: Optimal 4th order from R. Mclachlan, 4 drifts/4 kicks per step
%               6: Yoshida 6th order, 8 drifts/7 kicks per step
%
% [ROUT,LOST] = GPUPASS(...)
%   Additionally return information on lost particles
%
%   LOST        1x1 structure with the following fields:
%               turn                1xN vector, turn number where the particle is lost
%               element             1xN vector, element number where the particle is lost
%               coordinates_at_loss 6xN matrix, coordinates at the entrance of the
%                                   element where the particle was lost
%
% See also  RINGPASS, LINEPASS
error('at:missingMex','missing MEX file.');
end
