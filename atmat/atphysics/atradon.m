function [ring,radelems,cavities,energy]=atradon(ring1,varargin)
%ATRADON			switches RF and radiation on
%
%RING2=ATRADON(RING,CAVIPASS,BENDPASS,QUADPASS)
%   Changes passmethods to get RF cavity acceleration and radiation
%   damping. ATRADON also sets the "Energy" field on the modified elements,
%   looking for the machine energy in:
%       1) 1st 'RingParam' element
%       2) 1st 'RFCavity' element
%       3) field "E0" of the global variable "GLOBVAL"
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
%[RING2,RADINDEX,CAVINDEX,ENERGY]=ATRADON(...)
%           returns the index of radiative elements and cavities + energy

[cavipass,bendpass,quadpass]=parseargs({'CavityPass','auto',''},varargin);

ring=ring1;

energy=atenergy(ring);
if ~isempty(cavipass)
    cavities=atgetcells(ring,'Frequency');
    if ~any(cavities)
        warning('AT:atradon:NoCavity', 'No cavity found in the structure');
    end
    ring(cavities)=changepass(ring(cavities),cavipass,energy);
else
    cavities=false(size(ring));
end

if ~isempty(bendpass)
    isdipole=@(elem,bangle) bangle~=0;
    dipoles=atgetcells(ring,'BendingAngle',isdipole);
    if sum(dipoles) <= 0
        warning('AT:atradon:NoBend', 'No dipole in the structure');
    end
    ring(dipoles)=changepass(ring(dipoles),bendpass,energy);
else
    dipoles=false(size(ring));
end

if ~isempty(quadpass)
    isquadrupole=@(elem,polyb) length(polyb) >= 2 && polyb(2)~=0;
    quadrupoles=atgetcells(ring,'PolynomB',isquadrupole) & ~dipoles;
    if sum(quadrupoles) <= 0
        warning('AT:atradon:NoQuad', 'No quadrupole in the structure');
    end
    ring(quadrupoles)=changepass(ring(quadrupoles),quadpass,energy);
else
    quadrupoles=false(size(ring));
end

radelems=dipoles|quadrupoles;

disp(['Cavities located at position ' num2str(find(cavities)')]);
disp([num2str(sum(radelems)) ' elements switched to include radiation']);

    function newline=changepass(line,newpass,nrj)
    if strcmp(newpass,'auto')
        passlist=strrep(atgetfieldvalues(line,'PassMethod'),'Pass','RadPass');
    else
        passlist=repmat({newpass},size(line));
    end
    newline=cellfun(@newelem,line,passlist,'UniformOutput',false);

        function elem=newelem(elem,newpass)
            elem.PassMethod=newpass;
            elem.Energy=nrj;
        end
    end

end
