function energy=atenergy(ring)
%ENERGY=ATENERGY(RING) Gets the RING energy
%   ATENERGY looks for the machine energy in:
%       1) the 1st 'RingParam' element
%       2) the 1st 'RFCavity' element
%       3) the field "E0" of the global variable "GLOBVAL"
%

global GLOBVAL

params=atgetcells(ring,'Class','RingParam');
cavities=atgetcells(ring,'Frequency');
if any(params)
    energy=ring{find(params,1)}.Energy;
elseif any(cavities) && isfield(ring{find(cavities,1)},'Energy')
    energy=ring{find(cavities,1)}.Energy;
elseif isfield(GLOBVAL,'E0')
    energy=GLOBVAL.E0;
else
    error('AT:NoEnergy',...
        'Energy not defined (searched in ''RingParam'',''RFCavity'',GLOBVAL.E0)');
end
end

