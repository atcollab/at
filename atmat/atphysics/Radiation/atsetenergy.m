function newring = atsetenergy(ring,E)
% ATSETENERGY sets the Energy field in all elements.
% If no such field exists, it creates it.
%
% newring = atsetenergy(ring,Energy)
%
%   ring: an AT ring.
%   Energy: Value to set the Energy field. Units: eV
%   newring: new AT ring with Energy field set.
%
% Example:
%   Set the energy of the elements in RING to 3 GeV.
%   NEWRING = atsetenergy(RING,3e9)
%
% See also atenergy
%
newring=ring;
for j=1:length(ring)
    newring{j}.Energy = E;
end
