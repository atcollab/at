function [ ringID ] = add_ID( ring, IDelem )
%add_ID adds IDelem at the beginning and subtracts the length from
%   the surrounding straights
 
ringID={IDelem,ring{:}};
IDlen=IDelem.Length;

driftInd=findcells(ringID,'FamName','SDHI');
ringID{driftInd(1)}.Length=ringID{driftInd(1)}.Length-IDlen/2;
ringID{driftInd(end)}.Length=ringID{driftInd(end)}.Length-IDlen/2;

end

