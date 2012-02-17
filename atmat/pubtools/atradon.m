function [ring2,radindex,cavindex]=atradon(ring,cavipass,bendpass,quadpass)
%ATRADON			switches RF and radiation on
%
%RING2=ATRADON(RING,CAVIPASS,BENDPASS,QUADPASS)
%
%RING:		initial AT structure
%CAVIPASS:	pass method for cavities (default ThinCavityPass)
%BENDPASS:	pass method for cavities (default BndMPoleSymplectic4RadPass)
%QUADPASS:	pass method for cavities (default: nochange)
%
%[RING2,RADINDEX,CAVINDEX]=ATRADON(...) returns the index of radiative elements
%					 and cavities

global GLOBVAL
if nargin < 4, quadpass=[]; end
if nargin < 3, bendpass=[]; end
if nargin < 2
cavipass='CavityPass';
bendpass='BndMPoleSymplectic4RadPass';
end

if ~isfield(GLOBVAL,'E0')
   error('Energy not defined in GLOBVAL');
end

if ~isempty(cavipass)
   cavindex=findcells(ring,'Frequency');
   if isempty(cavindex)
      warning('No cavity found in the structure')
   else
      disp(['Cavities located at position ' num2str(cavindex)]);
   end
   ring2=setcellstruct(ring,'PassMethod',cavindex,cavipass);
else
   ring2=ring;
end

if ~isempty(bendpass)
   dipindex=findcells(ring2,'BetaCode','DI');
   if isempty(dipindex), warning('No dipole found in the structure'); end
   ring2=setcellstruct(ring2,'PassMethod',dipindex,bendpass);
else
   dipindex=[];
end

if ~isempty(quadpass)
   quadindex=findcells(ring2,'BetaCode','QP');
   if isempty(dipindex), warning('No quadrupole found in the structure'); end
   ring2=setcellstruct(ring2,'PassMethod',quadindex,quadpass);
else
   quadindex=[];
end

radindex=sort([dipindex quadindex]);
disp([num2str(length(radindex)) ' elements switched to include radiation']);
