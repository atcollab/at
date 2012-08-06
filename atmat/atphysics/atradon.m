function [ring2,radindex,cavindex]=atradon(ring,varargin)
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
quadpass='';
bendpass='BndMPoleSymplectic4RadPass';
cavipass='CavityPass';
if nargin >= 4, quadpass=varargin{3}; end
if nargin >= 3, bendpass=varargin{2}; end
if nargin >= 2, cavipass=varargin{1}; end

if ~isfield(GLOBVAL,'E0')
    error('Energy not defined in GLOBVAL');
end

if ~isempty(cavipass)
    cavities=atgetcells(ring,'Frequency');
    cavindex=find(cavities)';
    if sum(cavities) <= 0
        warning('AT:atradon:NoCavity',...
            'No cavity found in the structure')
    else
        disp(['Cavities located at position ' num2str(cavindex)]);
    end
    ring2=setcellstruct(ring,'PassMethod',cavities,cavipass);
else
    ring2=ring;
end

if ~isempty(bendpass)
    isdipole=@(elem,field) elem.(field)~=0;
    dipoles=atgetcells(ring2,'BendingAngle',isdipole);
    if sum(dipoles) <= 0, warning('AT:atradon:NoBend',...
            'No dipole found in the structure'); end
    ring2=setcellstruct(ring2,'PassMethod',dipoles,bendpass);
else
    dipoles=false(size(ring2));
end

if ~isempty(quadpass)
    isquadrupole=@(elem,field) length(elem.(field)) >= 2 && elem.(field)(2)~=0;
    quadrupoles=atgetcells(ring2,'PolynomB',isquadrupole) && ~dipoles;
    if sum(quadrupoles) <= 0, warning('AT:atradon:NoQuad',...
            'No quadrupole found in the structure'); end
    ring2=setcellstruct(ring2,'PassMethod',quadrupoles,quadpass);
else
    quadrupoles=false(size(ring2));
end

radindex=find(dipoles | quadrupoles)';
disp([num2str(length(radindex)) ' elements switched to include radiation']);
