function atclass=getclass_6d(elem)
%GETCLASS_6D    Private. Guess class for 6d motion
%
% Returns the element class even without a Class field

if isfield(elem,'Class')
    atclass=elem.Class;
elseif isfield(elem,'BendingAngle')
    atclass='Bend';
elseif isfield(elem,'PolynomB') && elem.Length > 0
    loworder=find(abs(elem.PolynomB(2:end))~=0,1);
    if isempty(loworder)
        atclass='Drift';
    elseif loworder==1
        atclass='Quadrupole';
    elseif loworder==2
        atclass='Sextupole';
    elseif loworder==3
        atclass='Octupole';
    else
        atclass='Multipole';
    end
elseif isfield(elem,'Frequency')
    atclass='RFCavity';
elseif isfield(elem,'Lmatp')
    atclass='QuantDiff';
else
    atclass='Other';
end

end
