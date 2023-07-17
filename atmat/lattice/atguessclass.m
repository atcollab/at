function atclass = atguessclass(elem, varargin)
%ATGUESSCLASS Tries to determine the class of an element
%ATCLASS=ATGUESSCLASS(ATELEM) Tries to determine the class of an element
%
%  INPUTS
%  1. elem - AT element
%
%
%  NOTES
%  1. atclass=atguessclass(atelem,'useclass')
%  By default, ATGUESSCLASS will default "useless" elements (PolynopmB==0)
%  to 'Drift' or 'Marker', depending on 'Length'. When specifying
%  'UseClass', it it will preserve the 'Class' field for those elements.
%
%  See also atwritem

useclass=any(strcmpi(varargin,'UseClass'));

if isfield(elem,'BendingAngle')
    atclass='Bend';
elseif isfield(elem,'Frequency')
    atclass='RFCavity';
elseif isfield(elem,'KickAngle')
    atclass='Corrector';
elseif isfield(elem,'Periodicity')
    atclass='RingParam';
elseif isfield(elem,'Limits')
    atclass='Aperture';
elseif isfield(elem,'PolynomB')
    if useclass && isfield(elem,'Class')
        atclass=elem.Class;
    elseif isfield(elem,'Length') && elem.Length~=0
        maxorder=elem.MaxOrder+1;
        loworder=find(abs(elem.PolynomB(2:maxorder))~=0,1);
        if isempty(loworder)
            atclass='Drift';
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
elseif isfield(elem,'Lmatp')
    atclass='QuantDiff';
elseif isfield(elem,'M66')
    if isfield(elem,'Tijk')
        atclass='MatrixTijkPass';
    else
        atclass='Matrix66';
    end
elseif isfield(elem,'xtable')
    atclass='InsertionDeviceKickMap';
elseif isfield(elem,'Length') && elem.Length~=0
    atclass='Drift';
else
    if useclass && isfield(elem,'Class')
        atclass=elem.Class;
    else
        atclass='Marker';
    end
end

end
