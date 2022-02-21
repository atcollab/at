function [energy,nbper,voltage,harmnumber,U0]=atenergy(ring)
%ATENERGY Gets the lattice energy
%
%  ENERGY=ATENERGY(RING)
%  [ENERGY,PERIODS]=atenergy(RING)
%  [ENERGY,PERIODS,VOLTAGE,HARMNUMBER]=atenergy(RING)
%  [ENERGY,PERIODS,VOLTAGE,HARMNUMBER,U0]=atenergy(RING)
%
% Warning: To get ENERGY, PERIODS and HARMNUMBER, use atGetRingProperties
%          To get U0, use atgetU0
%
%   RING        Ring structure
%
%   ENERGY      Ring energy
%       ATENERGY looks for the machine energy in:
%           1) the 1st 'RingParam' element
%           2) the 'RFCavity' with the lowest frequency
%           3) the field "E0" of the global variable "GLOBVAL"
%           4) The field "Energy" in any element
%   PERIODS     Number of periods
%   VOLTAGE     Total RF voltage for the main cavities. The main cavities
%               are the ones with the lowest frequency
%   HARMNUMBER  Harmonic number. Computed from the frequency of the main cavities
%   U0          Total energy loss per turn
%
%  See also atGetRingProperties atgetU0 atsetcavity

global GLOBVAL %#ok<GVMIS> 
s=warning;                          % Save the warning state
warning('Off','AT:NoRingParam');    % Disable warning
props = atGetRingProperties(ring);  % Get ring propeties
warning(s);                         % Restore the warning state
energy = props.Energy;
if ~isfinite(energy)
    if isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    else
        error('AT:NoEnergy',...
                ['Energy not defined (searched in ''RingParam'',''RFCavity'', ',...
                'the ''Energy'' field of each element, GLOBVAL.E0)\n' ...
                'You can set it with:\n>> ring=atSetRingProperties(ring,''energy'',energy)\n'...
                'or by using the global variable GLOBVAL:\n>> GLOBVAL.E0=energy;'...
                ]);
    end
end
nbper = props.Periodicity;

if nargout >= 3
    if isfield(props, 'cavpts')
        maincavs = ring(props.cavpts);
    else
        maincavs = findmaincavs(ring);
    end
    if ~isempty(maincavs)
        voltage=nbper*sum(atgetfieldvalues(maincavs,'Voltage'));
    elseif nargout >= 5
        voltage=NaN;
    else
        error('AT:NoCavity','No cavity element in the ring');
    end
end
if nargout >= 4
    if isfield(props,'HarmNumber')
        harmnumber = props.HarmNumber;
    else
        harmnumber=NaN;
    end
end
if nargout >= 5
    U0 = atgetU0(ring);
end

    function cavs=findmaincavs(ring)
        cavities = atgetcells(ring,'Frequency');
        if any(cavities)
            freqs = atgetfieldvalues(ring(cavities), 'Frequency');
            [~,~,ic]=unique(freqs);
            cavities(cavities) = (ic == 1);
            cavs = ring(cavities);
        else
            cavs={};
        end
    end

end
