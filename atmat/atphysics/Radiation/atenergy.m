function [energy,nbper,voltage,harmnumber,U0]=atenergy(ring)
%ATENERGY Gets the lattice energy
%
%  ENERGY=ATENERGY(RING) Gets the RING energy
%   ATENERGY looks for the machine energy in:
%       1) the 1st 'RingParam' element
%       2) the 1st 'RFCavity' element
%       3) the field "E0" of the global variable "GLOBVAL"
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
%    1. Cavities in the ring must have the Class RFCavity.
%    2. Check for 2 pi phase advance if more than 2 ouput arguments
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
dipoles  = atgetcells(ring(:,1),'BendingAngle');
theta    = atgetfieldvalues(ring(dipoles),'BendingAngle');

if any(params) % case RingParam is defined in the lattice
    parmelem=ring{find(params,1)};
    energy=parmelem.Energy;
    if nargout >= 2
        nbper=parmelem.Periodicity;
    end
else % else look for energy in cavity or GLOBVAL
    if any(cavities) && isfield(ring{find(cavities,1)},'Energy')
        energy=ring{find(cavities,1)}.Energy;
    elseif isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    else
        error('AT:NoEnergy',...
            'Energy not defined (searched in ''RingParam'',''RFCavity'',GLOBVAL.E0)');
    end
    if nargout >= 2 % get number of periods
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
        harmnumber=nbper*atgetfieldvalues(ring(find(cavities,1)),'HarmNumber');
    elseif nargout >= 5
        voltage=NaN;
        harmnumber=NaN;
    else
        error('AT:NoCavity','No cavity element in the ring');
    end
end

% get energy loss per turn
if nargout >= 5 
    % Losses = Cgamma/2/pi*EGeV^4*I2
    cgamma=4e9*pi*PhysConstant.classical_electron_radius.value/3/...
        PhysConstant.electron_mass_energy_equivalent_in_MeV.value^3; % [m/GeV^3]
    lendp=atgetfieldvalues(ring(dipoles),'Length');
    losses=atgetfieldvalues(ring(atgetcells(ring,'I2')),'I2');
    I2=nbper*(sum(abs(theta.*theta./lendp))+sum(losses));            % [m-1]
    U0=cgamma/2/pi*(energy*1.e-9)^4*I2*1e9;                          % [eV]
end
end
