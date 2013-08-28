function [ring2,radelems,cavities,energy]=atradon(ring,varargin)
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

global GLOBVAL
[cavipass,bendpass,quadpass]=decodeargs({'CavityPass','auto',''},varargin);

params=atgetcells(ring,'Class','RingParam');
cavities=atgetcells(ring,'Frequency');
if any(params)
    energy=ring{find(params,1)}.Energy;
elseif any(cavities) && isfield(ring{find(cavities,1)},'Energy')
    energy=ring{find(cavities,1)}.Energy;
elseif isfield(GLOBVAL,'E0')
    energy=GLOBVAL.E0;
else
    error('AT:atradon:NoGLOBVAL',...
        'Energy not defined (searched in ''RingParam'',''RFCavity'',GLOBVAL.E0)');
end

ring2=ring;
if ~isempty(cavipass)
    if ~any(cavities)
        warning('AT:atradon:NoCavity', 'No cavity found in the structure');
    end
    ring2(cavities)=changepass(ring2(cavities),cavipass);
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

radelems=dipoles|quadrupoles;

disp(['Cavities located at position ' num2str(find(cavities)')]);
disp([num2str(sum(radelems)) ' elements switched to include radiation']);

    function newline=changepass(line,newpass)
    if strcmp(newpass,'auto')
        passlist=strrep(atgetfieldvalues(line,'PassMethod'),'Pass','RadPass');
    else
        passlist=repmat({newpass},size(line));
    end
    newline=cellfun(@newelem,line,passlist,'UniformOutput',false);

        function elem=newelem(elem,newpass)
            elem.PassMethod=newpass;
            elem.Energy=energy;
        end
    end

end
