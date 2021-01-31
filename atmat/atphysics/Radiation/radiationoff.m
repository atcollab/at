%RADIATIONOFF turns classical radiation  OFF
%  Switch all magnets currently set to use pass-methods
%  'BndMPoleSymplectic4RadPass' and  'StrMPoleSymplectic4RadPass'
%  to their equivalents without radiation
%  'BndMPoleSymplectic4Pass' and  'StrMPoleSymplectic4Pass'
%	
%  NOTES:
%    1. Deprecated function, use atradoff instead
%
%   See also RADIATIONON, CAVITYON, CAVITYOFF, ATRADON, ATRADOFF


if ~evalin('base','exist(''THERING'')') || ~evalin('base','~isempty(whos(''global'',''THERING''))')
   error('Global variable THERING could not be found');
end
localindex = findcells(THERING,'PassMethod','StrMPoleSymplectic4RadPass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'StrMPoleSymplectic4Pass');
totalswitched = length(localindex);

localindex = findcells(THERING,'PassMethod','BndMPoleSymplectic4RadPass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'BndMPoleSymplectic4Pass');
totalswitched = totalswitched + length(localindex);

localindex = findcells(THERING,'PassMethod','BndMPoleSymplectic4E2RadPass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'BndMPoleSymplectic4E2Pass');
totalswitched = totalswitched + length(localindex);

disp(['PassMethod was changed to NOT include radiation in ',num2str(totalswitched),  ' elements'])     
clear localindex totalswitched
