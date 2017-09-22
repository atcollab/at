function fnames = entrancefields()
%ENTRANCEFIELDS() Return the list of field names affecting the element entrance

persistent fnms
if isempty(fnms)
    fnms={'T1','R1','EntranceAngle','FringeInt1',...
        'FringeBendEntrance','FringeQuadEntrance'};
end
fnames=fnms;
end

