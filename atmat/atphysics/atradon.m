function [ring,radelemIndex,cavitiesIndex,energy]=atradon(ring1,varargin)
%ATRADON switches RF and radiation on
%
%  [RING2,RADINDEX,CAVINDEX,ENERGY] = ATRADON(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to get RF cavity acceleration and radiation
%    damping. ATRADON also sets the "Energy" field on the modified elements,
%    looking for the machine energy in:
%       1) 1st 'RingParam' element
%       2) 1st 'RFCavity' element
%       3) field "E0" of the global variable "GLOBVAL"
%
%  INPUTS
%  1. RING	     initial AT structure
%  2. CAVIPASS   pass method for cavities (default ThinCavityPass)
%                '' makes no change,
%  3. BENDPASS   pass method for bending magnets. Special values:
%                '' makes no change,
%                'auto' will substitute 'Pass' with 'RadPass' in any method
%                (default: 'auto')
%  4. QUADPASS   pass method for quadrupoles
%                '' makes no change,
%                'auto' will substitute 'Pass' with 'RadPass' in any method
%                (default: '')
%
%  OUPUTS
%  1. RING2     Output ring
%  2. RADINDEX  Indices of elements with radiation
%  3. CAVINDEX  Indices of cavities
%
%  See also ATRADOFF, ATCAVITYON, ATCAVITYOFF


[cavipass,bendpass,quadpass]=parseargs({'CavityPass','auto',''},varargin);

ring=ring1;

energy=atenergy(ring);
if ~isempty(cavipass)
    cavitiesIndex=atgetcells(ring,'Frequency');
    if any(cavitiesIndex)
        ring(cavitiesIndex)=changepass(ring(cavitiesIndex),cavipass,energy);
    end
else
    cavitiesIndex=false(size(ring));
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

radelemIndex=dipoles|quadrupoles;

if any(cavitiesIndex)
    atdisplay(1,['Cavities located at position ' num2str(find(cavitiesIndex)')]);
else
    atdisplay(1,'No cavity');
end
atdisplay(1,[num2str(sum(radelemIndex)) ' elements switched to include radiation']);

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
