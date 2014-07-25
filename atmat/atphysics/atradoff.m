function [ring,radelems]=atradoff(ring1,varargin)
%ATRADOFF			switches radiation off
%
%RING2=ATRADOFF(RING,CAVIPASS,BENDPASS,QUADPASS)
%   Changes passmethods to turn off radiation
%   damping. 
%
%RING:		initial AT structure
%BENDPASS:	pass method for bending magnets. Special values:
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: 'auto')
%QUADPASS:	pass method for quadrupoles
%           '' makes no change,
%           'auto' wille substitute 'Pass' with 'RadPass' in any method
%           (default: '')
%
%[RING2,radelems,cavities]=ATRADOFF(...)
%           returns the index of radiative elements and cavities

[bendpass,quadpass]=parseargs({'auto',''},varargin);

ring=ring1;

energy=atenergy(ring);


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

disp([num2str(sum(radelems)) ' elements with radiation switched off']);

    function newline=changepass(line,newpass,nrj)
    if strcmp(newpass,'auto')
        passlist=strrep(atgetfieldvalues(line,'PassMethod'),'RadPass','Pass');
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
