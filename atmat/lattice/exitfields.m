function fnames = exitfields(varargin)
%EXITFIELDS() Return the list of field names affecting the element exit
%
% Optional arguments:
% 'KeepAxis', if present, rotations translations are excluded from list 
%
%see also: entrancefields atdivelem

rottrasl=getflag(varargin,'KeepAxis');

start=1;

if rottrasl
   start=3; % remove R1 T1 from list of field at entrance
end

persistent fnms
if isempty(fnms)
    fnms={'T2','R2','ExitAngle','FringeInt2',...
        'FringeBendExit','FringeQuadExit'};
end
fnames=fnms(start:end);
end

