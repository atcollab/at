%RADIATIONON turns classical radiation  ON
% Switch all magnets currently set to use pass-methods
% 'BndMPoleSymplectic4Pass' and  'StrMPoleSymplectic4Pass'
% to their equivalents with classical radiation
% 'BndMPoleSymplectic4RadPass' and  'StrMPoleSymplectic4RadPass'
%	
%   See also RADIATIONOFF, CAVITYON, CAVITYOFF


if ~isglobal(THERING)
   error('Global variable THERING could not be found');
end
localindex = findcells(THERING,'PassMethod','StrMPoleSymplectic4Pass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'StrMPoleSymplectic4RadPass');
totalswitched = length(localindex);

localindex = findcells(THERING,'PassMethod','BndMPoleSymplectic4Pass');
THERING = setcellstruct(THERING,'PassMethod',localindex, 'BndMPoleSymplectic4RadPass');
totalswitched = totalswitched + length(localindex);

disp(['PassMethod was changed to include radiation in ',num2str(totalswitched),  ' elements']) 
clear localindex