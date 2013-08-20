function [ring2,radindex,cavindex]=atradon(ring,varargin)
%ATRADON			switches RF and radiation on
%
%RING2=ATRADON(RING,CAVIPASS,BENDPASS,QUADPASS)
%
%RING:		initial AT structure
%CAVIPASS:	pass method for cavities (default ThinCavityPass)
%BENDPASS:	pass method for bending magnets. Special values:
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: 'auto')
%QUADPASS:	pass method for quadrupoles
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: '')
%
%[RING2,RADINDEX,CAVINDEX]=ATRADON(...) returns the index of radiative elements
%					 and cavities

global GLOBVAL
args={'CavityPass','auto',''};
args(1:length(varargin))=varargin;
[cavipass,bendpass,quadpass]=deal(args{:});

if ~isfield(GLOBVAL,'E0')
    error('AT:atradon:NoGLOBVAL','Energy not defined in GLOBVAL');
end

ring2=ring;
if ~isempty(cavipass)
    cavities=atgetcells(ring,'Frequency');
    cavindex=find(cavities)';
    if sum(cavities) <= 0
        warning('AT:atradon:NoCavity', 'No cavity found in the structure');
    else
        disp(['Cavities located at position ' num2str(cavindex)]);
    end
    ring2(cavities)=atsetfieldvalues(ring(cavities),'PassMethod',cavipass);
else
    ring2=ring;
end

if ~isempty(bendpass)
    isdipole=@(elem,bangle) bangle~=0;
    dipoles=atgetcells(ring2,'BendingAngle',isdipole);
    if sum(dipoles) <= 0
        warning('AT:atradon:NoBend', 'No dipole in the structure');
    end
    ring2(dipoles)=changepass(ring2(dipoles),bendpass);
else
    dipoles=false(size(ring2));
end

if ~isempty(quadpass)
    isquadrupole=@(elem,polyb) length(polyb) >= 2 && polyb(2)~=0;
    quadrupoles=atgetcells(ring2,'PolynomB',isquadrupole) & ~dipoles;
    if sum(quadrupoles) <= 0
        warning('AT:atradon:NoQuad', 'No quadrupole in the structure');
    end
    ring2(quadrupoles)=changepass(ring2(quadrupoles),quadpass);
else
    quadrupoles=false(size(ring2));
end

radindex=find(dipoles | quadrupoles)';
disp([num2str(length(radindex)) ' elements switched to include radiation']);

    function newline=changepass(line,newpass)
    if strcmp(newpass,'auto')
        newpass=strrep(atgetfieldvalues(line,'PassMethod'),'Pass','RadPass');
    end
    newline=atsetfieldvalues(line,'PassMethod',newpass);
    end

end
