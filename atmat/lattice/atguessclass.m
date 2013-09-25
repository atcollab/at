function atclass = atguessclass(elem)
%ATCLASS=ATGUESSCLASS(ATELEM) Tries to determine the class of an element
%
if isfield(elem,'BendingAngle')
    atclass='Bend';
elseif isfield(elem,'Frequency')
    atclass='RFCavity';
elseif isfield(elem,'KickAngle')
    atclass='Corrector';
elseif isfield(elem,'PolynomB')
    if isfield(elem,'Length') && elem.Length~=0
        loworder=find(abs(elem.PolynomB(2:end))~=0,1);
        if loworder==1
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
    if isfield(elem,'Length') && elem.Length~=0
        atclass='Drift';
    else
        atclass='Marker';
    end
end
end
