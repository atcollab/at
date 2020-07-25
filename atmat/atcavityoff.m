function [ring, ATCavityIndex]=atcavityoff(ring)
%ATCAVITYOFF switches RF cavities off
%
%  [RING2, CAVITIESINDEX] = ATCAVITYOFF(RING)
%    Changes passmethods to turn off radiation
%   damping.
%
%  INPUTS:
%  1. RING      initial AT structure
%
%  OUPUTS
%  1. RING2          output ring with cavities switched off
%  2. CAVITIESINDEX  indices cavities
%
%  See also ATCAVITYON, ATRADON, ATRADOFF

%
%% Written by Laurent S. Nadolski

% Look for cavities
ATCavityIndex = findcells(ring, 'Frequency');

% Return if no cavity found
if isempty(ATCavityIndex)
    atdisplay(1,'AT:atcavityoff: No cavities were found in the lattice.');
    return
end

% Make column vector 
ATCavityIndex =ATCavityIndex(:)';

% Loops over cavity elements
for iCavity =ATCavityIndex        
    % Based on cavity length, decide of the passmethod to switch off the cavities
    if ring{iCavity}.Length == 0
        ring{iCavity}.PassMethod = 'IdentityPass';
    else
        ring{iCavity}.PassMethod = 'DriftPass';
    end
end

end
