function [ring, ATCavityIndex]=atcavityon(ring,cavityPass)
%ATRADON switches RF cavities on
%
%  [RING2,CAVINDEX]=ATCAVITYON(RING,CAVITYPASS)
%    Changes passmethods to get RF cavity acceleration and radiation
%    damping. 
%
%  INPUTS
%  1. RING	     	initial AT structure
%  2. CAVITYPASS	customed passmethod for cavities (default CavityPass)      
%
%  OUPUTS
%  1. RING2          output ring with cavities off
%  2. CAVITIESINDEX  indices of radiative elements and cavities
%  3. ENERGY         energy
%
%  See also ATRADOFF, ATRADON, ATCAVITYOFF

%
%% Written by Laurent S. Nadolski

if nargin == 1 
    cavityPass = 'On';
end

% Look for cavity location
ATCavityIndex = findcells(ring, 'Frequency');

% Return if no cavity found
if isempty(ATCavityIndex)
    atdisplay(1,'AT:atcavityon: No cavities were found in the lattice.');
    return
end

% Make column vector 
ATCavityIndex =ATCavityIndex(:)';

% Loops over cavity elements
for iCavity =ATCavityIndex
    
    if size(cavityPass,1) == 1
        CavityString = deblank(cavityPass);
    elseif size(cavityPass,1) == length(ATCavityIndex)
        CavityString = deblank(cavityPass(iCavity,:));
    else
        error('AT:atcavityoff: Number of rows in the input string must be 1 row or equal to the number of cavities.');
    end
    
    % Set cavity pass method    
    if strcmpi(CavityString,'On')
        ring{iCavity}.PassMethod = 'CavityPass';        
    else % Custom passmethod
        ring{iCavity}.PassMethod = CavityString;
    end
end

end
