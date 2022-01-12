function [energy,nbper,voltage,harmnumber,U0]=atenergy(ring)
%ATENERGY Gets the lattice energy
%
%  ENERGY=ATENERGY(RING) Gets the RING energy
%   ATENERGY looks for the machine energy in:
%       1) the 1st 'RingParam' element
%       2) the 1st 'RFCavity' element
%       3) the field "E0" of the global variable "GLOBVAL"
%       4) The field "Energy" in any element
%  
%  INPUTS
%    1. ring  Ring structure
%
%  OUPUTS
%    1. energy     Ring energy in eV
%    2. nbper      Number of periods to make 2pi phase advance
%    3. voltage    Total RF voltage  
%    4. harmnumber Harmonic number
%    5. U0         Energy loss per turn in eV
%
%  NOTES
%    1. voltage and harmnumber assume a single RF system: fundamental
%       frequency, no harmonic system. More prcisely:
%
%       voltage is the sum of the voltages of all cavities,
%       harmnuber is computed with the frequency of the 1st cavity in the ring
%
%  EXAMPLES
%    1. [ENERGY,PERIODS]=atenergy(ring) also outputs the number of periods
%    2. [ENERGY,PERIODS,VOLTAGE,HARMNUMBER]=atenergy(ring) also outputs the harmonic number
%    3. [ENERGY,PERIODS,VOLTAGE,HARMNUMBER,U0]=atenergy(ring) also outputs total losses in eV
%
%  See also atgetU0 atsetcavity atenergy

global GLOBVAL
TWO_PI_ERROR = 1.e-4;
 
params   = atgetcells(ring(:,1),'Class','RingParam');
cavities = atgetcells(ring(:,1),'Frequency');
        
if any(params) % case RingParam is defined in the lattice
    parmelem=ring{find(params,1)};
    energy=parmelem.Energy;
    if nargout >= 2
        nbper=parmelem.Periodicity;
    end
else % else look for energy in cavity or GLOBVAL
    E0s = atgetfieldvalues(ring,'Energy');
    E0s = E0s(isfinite(E0s));         % Discard undefined values
    if any(cavities) && isfield(ring{find(cavities,1)},'Energy')
        energy=ring{find(cavities,1)}.Energy;
    elseif isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    elseif length(unique(E0s)) == 1
        energy = unique(E0s);
    elseif length(unique(E0s)) > 1
        error('AT:NoEnergy','Energy field not equal for all elements')
    else
        error('AT:NoEnergy',...
            ['Energy not defined (searched in '...
            '''RingParam'',''RFCavity'',GLOBVAL.E0,',...
            ' field ''Energy'' of each element)']);
    end
    if nargout >= 2 % get number of periods
        dipoles  = atgetcells(ring(:,1),'BendingAngle');
        theta    = atgetfieldvalues(ring(dipoles),'BendingAngle');
        if size(ring,2) > 1
            nbper=size(ring,2);
        else % Guess number of periods based on the total bending angle
            nbp=2*pi/sum(theta);
            nbper=round(nbp);
            if ~isfinite(nbp)
                warning('AT:WrongNumberOfCells','No bending in the cell, ncells set to 1');
                nbper=1;
            elseif abs(nbp-nbper) > TWO_PI_ERROR
                warning('AT:WrongNumberOfCells','non integer number of cells: ncells = %g -> %g',nbp,nbper);
            end
        end
    end
end

% Get total cavity voltage and harmonic number
if nargout >= 3 
    if any(cavities) % sum over cavity number
        voltage=nbper*sum(atgetfieldvalues(ring(cavities),'Voltage'));
        %harmnumber=nbper*atgetfieldvalues(ring(find(cavities,1)),'HarmNumber');
        Frf=atgetfieldvalues(ring(find(cavities,1)),'Frequency');
        L0 = nbper*findspos(ring,length(ring)+1); % design length [m]
        C0 = PhysConstant.speed_of_light_in_vacuum.value; % speed of light [m/s]
        harmnumber = round(Frf*L0/C0);

    elseif nargout >= 5
        voltage=NaN;
        harmnumber=NaN;
    else
        error('AT:NoCavity','No cavity element in the ring');
    end
end

% get energy loss per turn
if nargout >= 5
    U0 = atgetU0(ring);
end

end
