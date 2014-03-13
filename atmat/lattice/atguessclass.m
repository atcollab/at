function atclass = atguessclass(elem, varargin)
%ATCLASS=ATGUESSCLASS(ATELEM) Tries to determine the class of an element
%   ELEM    AT element
%
%ATCLASS=ATGUESSCLASS(ATELEM,'UseClass')
%   By default, ATGUESSCLASS will default "useless" elements (PolynopmB==0)
%   to 'Drift' or 'Marker', depending on 'Length'. When specifying
%   'UseClass', it it will preserve the 'Class' field for those elements.
%

useclass=any(strcmpi(varargin,'UseClass'));

if isfield(elem,'BendingAngle')
    atclass='Bend';
elseif isfield(elem,'Frequency')
    atclass='RFCavity';
elseif isfield(elem,'KickAngle')
    atclass='Corrector';
elseif isfield(elem,'PolynomB')
    if isfield(elem,'Length') && elem.Length~=0
        loworder=find(abs(elem.PolynomB(2:end))~=0,1);
        if isempty(loworder)
            if useclass && isfield(elem,'Class')
                atclass=elem.Class;
            else
                atclass=defaultclass(elem);
            end
        elseif loworder==1
            atclass='Quadrupole';
        elseif loworder==2
            atclass='Sextupole';
        else
            atclass='Multipole';
        end
    else
        atclass='ThinMultipole';
    end
elseif isfield(elem,'Periodicity')
    atclass='RingParam';
else
    atclass=defaultclass(elem);
end

    function atc=defaultclass(elem)
    if isfield(elem,'Length') && elem.Length~=0
        atc='Drift';
    else
        atc='Marker';
    end
    end

end
