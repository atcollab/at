function [ring, ATCavityIndex]=atcavityon(ring,cavityPass)
%ATRADON switches RF cavities on
%
%  [RING2,CAVINDEX]=ATCAVITYON(RING,CAVITYPASS)
%    Changes cavity passmethods to get RF acceleration
%
%  INPUTS
%  1. RING	     	initial AT structure
%  2. CAVITYPASS	customed passmethod for cavities (default RFCavityPass)
%
%  OUPUTS
%  1. RING2          output ring with cavities off
%  2. CAVITIESINDEX  indices of cavities
%
%  See also ATCAVITYOFF, ATRADON, ATRADOFF

%% Written by Laurent S. Nadolski

if nargin <= 1
    cavityPass='RFCavityPass';
end

ATCavityIndex=atgetcells(ring,'Frequency');
ring=atradon(ring,cavityPass,'','','','','');

end
