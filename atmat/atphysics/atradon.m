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
%           '' makes no change,
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
    if any(cavities)
        ring(cavities)=changepass(ring(cavities),cavipass,energy);
    end
else
    cavities=false(size(ring));
end

if ~isempty(bendpass)
    isdipole=@(elem,bangle) bangle~=0;
    dipoles=atgetcells(ring,'BendingAngle',isdipole);
    if any(dipoles) > 0
        ring(dipoles)=changepass(ring(dipoles),bendpass,energy);
    end
else
    dipoles=false(size(ring));
end

if ~isempty(quadpass)
    isquadrupole=@(elem,polyb) length(polyb) >= 2 && polyb(2)~=0;
    quadrupoles=atgetcells(ring,'PolynomB',isquadrupole) & ~dipoles;
    if any(quadrupoles) > 0
        ring(quadrupoles)=changepass(ring(quadrupoles),quadpass,energy);
    end
else
    quadrupoles=false(size(ring));
end

radelems=dipoles|quadrupoles;

if any(cavities)
    atdisplay(1,['Cavities located at position ' num2str(find(cavities)')]);
else
    atdisplay(1,'No cavity');
end
atdisplay(1,[num2str(sum(radelems)) ' elements switched to include radiation']);

    function newline=changepass(line,newpass,nrj)
    if strcmp(newpass,'auto')
        passlist=atgetfieldvalues(line,'PassMethod');
        ok=cellfun(@(psm) isempty(strfind(psm,'RadPass')),passlist);
        passlist(ok)=strrep(passlist(ok),'Pass','RadPass');
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
