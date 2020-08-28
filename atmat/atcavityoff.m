function [ring, ATCavityIndex]=atcavityoff(ring)
%ATCAVITYOFF switches RF cavities off
%
%  [RING2, CAVITIESINDEX] = ATCAVITYOFF(RING)
%    Changes cavity passmethods to turn off acceleration
%
%  INPUTS:
%  1. RING      initial AT structure
%
%  OUPUTS
%  1. RING2          output ring with cavities switched off
%  2. CAVITIESINDEX  indices cavities
%
%  See also ATCAVITYON, ATRADON, ATRADOFF

%% Written by Laurent S. Nadolski

[ring,~,ATCavityIndex]=atradoff(ring,'auto','','','','');

end
