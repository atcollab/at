function [energy,nbper,voltage,harmnumber]=atenergy(ring)
%ENERGY=ATENERGY(RING) Gets the RING energy
%   ATENERGY looks for the machine energy in:
%       1) the 1st 'RingParam' element
%       2) the 1st 'RFCavity' element
%       3) the field "E0" of the global variable "GLOBVAL"
%
%[ENERGY,PERIODS]=ATENERGY(RING) also outputs the number of periods
%
%[ENERGY,PERIODS,VOLTAGE,HARMNUMBER]=ATENERGY(RING) also outputs the harmonic number

global GLOBVAL

params=atgetcells(ring(:,1),'Class','RingParam');
cavities=atgetcells(ring(:,1),'Frequency');
if any(params)
    parmelem=ring{find(params,1)};
    energy=parmelem.Energy;
    nbper=parmelem.Periodicity;
else
    if any(cavities) && isfield(ring{find(cavities,1)},'Energy')
        energy=ring{find(cavities,1)}.Energy;
    elseif isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    else
        error('AT:NoEnergy',...
        'Energy not defined (searched in ''RingParam'',''RFCavity'',GLOBVAL.E0)');
    end
    if size(ring,2) > 1
        nbper=size(ring,2);
    else
        nbp=2*pi/sum(atgetfieldvalues(ring(atgetcells(ring,'BendingAngle')),'BendingAngle'));
        nbper=round(nbp);
        if ~isfinite(nbp)
            warning('AT:WrongNumberOfCells','No bending in the cell, ncells set to 1');
            nbper=1;
        elseif abs(nbp-nbper) > 1.e-4
            warning('AT:WrongNumberOfCells','non integer number of cells: ncells = %g -> %g',nbp,nbper);
        end
    end
end

if nargout >= 3
    if any(cavities)
        voltage=nbper*sum(atgetfieldvalues(ring(cavities),'Voltage'));
        harmnumber=nbper*atgetfieldvalues(ring(find(cavities,1)),'HarmNumber');
    else
        error('AT:NoCavity','No cavity element in the ring');
    end
end
end
