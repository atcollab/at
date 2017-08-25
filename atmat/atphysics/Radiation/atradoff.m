function [ring2,radelemIndex,cavitiesIndex]=atradoff(ring1,varargin)
%ATRADOFF  switches radiation and cavity off
%
%   [RING2,RADINDEX,CAVINDEX] = ATRADOFF(RING,CAVIPASS,BENDPASS,QUADPASS)
%    Changes passmethods to turn off radiation damping. 
%
%  INPUTS
%  1. RING      initial AT structure
%  2. CAVIPASS  pass method for cavities (default IdentityPass)
%               '' makes no change,
%  3. BENDPASS  pass method for bending magnets. Special values:
%               '' makes no change,
%               'auto' wille substitute 'RadPass' with 'Pass' in any method
%               (default: 'auto')
%  4. QUADPASS  pass method for quadrupoles
%               '' makes no change,
%               'auto' wille substitute 'RadPass' with 'Pass' in any method
%               (default: '')
%
%   OUPUTS
%   1. RING2     Output ring
%   2. RADINDEX  Indices of elements with radiation
%   3. CAVINDEX  Indices of cavities
%
%  See also ATRADON, ATCAVITYON, ATCAVITYOFF

[cavipass,bendpass,quadpass]=parseargs({'IdentityPass','auto',''},varargin);

ring2=ring1;

if ~isempty(cavipass)
    cavitiesIndex=atgetcells(ring2,'Frequency');
    if ~any(cavitiesIndex)
        warning('AT:atradon:NoCavity', 'No cavity found in the structure');
    end
    ring2(cavitiesIndex)=changepass(ring2(cavitiesIndex),cavipass);
else
    cavitiesIndex=false(size(ring2));
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

radelemIndex=dipoles|quadrupoles;

disp(['Cavities located at position ' num2str(find(cavitiesIndex)')]);
disp([num2str(sum(radelemIndex)) ' elements with radiation switched off']);

    function newline=changepass(line,newpass)
    if strcmp(newpass,'auto')
        passlist=strrep(atgetfieldvalues(line,'PassMethod'),'RadPass','Pass');
    else
        passlist=repmat({newpass},size(line));
    end
    newline=cellfun(@newelem,line,passlist,'UniformOutput',false);

        function elem=newelem(elem,newpass)
            elem.PassMethod=newpass;
            %elem=rmfield(elem,'Energy');
        end
    end

end
