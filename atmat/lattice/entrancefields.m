function fnames = entrancefields(varargin)
%ENTRANCEFIELDS() Return the list of field names affecting the element entrance
%
% Optional arguments:
% 'KeepAxis', if present, rotations translations are excluded from list 
%
%see also: exitfields atdivelem

rottrasl=getflag(varargin,'KeepAxis');

start=1;

if rottrasl
   start=3; % remove R1 T1 from list of field at entrance
end

persistent fnms
if isempty(fnms)
    fnms={'T1','R1','EntranceAngle','FringeInt1',...
        'FringeBendEntrance','FringeQuadEntrance'};
end
fnames=fnms(start:end);
end

