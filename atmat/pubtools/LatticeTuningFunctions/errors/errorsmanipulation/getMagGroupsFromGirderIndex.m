function maggroups=getMagGroupsFromGirderIndex(r)
%
% output maggroups in r with indexes between GS and GE markers.
%
% maggroups is a cell array of magnet indexes describing a single magnet in
% reality, but sliced in the lattice
% a single magnet has the same MagNum value.
% 
%see also: UniformGirderErrors

indGS=find(atgetcells(r,'FamName','GS'));
indGE=find(atgetcells(r,'FamName','GE'));

maggroups=arrayfun(@(a,b)makegroup(a,b),indGS,indGE,'un',0);

return

function g=makegroup(a,b)
g=a:1:b;
return