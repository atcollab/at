function fnames = exitfields()
%EXITFIELDS() Return the list of field names affecting the element exit

persistent fnms
if isempty(fnms)
    fnms={'T2','R2','ExitAngle','FringeInt2',...
        'FringeBendExit','FringeQuadExit'};
end
fnames=fnms;
end

