%RADIATIONON turns classical radiation  ON
%
%  Switch all magnets currently set to use pass-methods
%  'BndMPoleSymplectic4Pass' and  'StrMPoleSymplectic4Pass'
%  to their equivalents with classical radiation
%  'BndMPoleSymplectic4RadPass' and  'StrMPoleSymplectic4RadPass'
%
%  NOTES:
%    1. Deprecated function, use atradon instead
%	
%   See also RADIATIONOFF, CAVITYON, CAVITYOFF, ATRADON, ATRADOFF


if ~evalin('base','exist(''THERING'')') || ~evalin('base','~isempty(whos(''global'',''THERING''))')
   error('Global variable THERING could not be found');
end

localindex = findcells(THERING,'PassMethod','StrMPoleSymplectic4Pass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'StrMPoleSymplectic4RadPass');
totalswitched = length(localindex);

localindex = findcells(THERING,'PassMethod','BndMPoleSymplectic4Pass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'BndMPoleSymplectic4RadPass');
totalswitched = totalswitched + length(localindex);

localindex = findcells(THERING,'PassMethod','BndMPoleSymplectic4E2Pass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'BndMPoleSymplectic4E2RadPass');
totalswitched = totalswitched + length(localindex);

disp(['PassMethod was changed to include radiation in ',num2str(totalswitched),  ' elements']) 
clear localindex