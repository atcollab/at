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
%    1. Check for 2 pi phase advance if more than 2 ouput arguments
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
E0s = atgetfieldvalues(ring,'Energy');
E0s = E0s(finite(E0s));         % Discard undefined values
        
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
    e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
    cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s
    e_radius=PhysConstant.classical_electron_radius.value;                  % m
    Cgamma=4.0E27*pi*e_radius/3/e_mass^3;                       % [m/GeV^3]
    Brho=sqrt(energy*energy - e_mass*e_mass)/cspeed;
    
    % Dipole radiation
    lendp=atgetfieldvalues(ring(dipoles),'Length');
    I2d=sum(abs(theta.*theta./lendp));
    % Wiggler radiation
    iswiggler=@(elem) strcmp(elem.Class,'Wiggler') && ~strcmp(elem.PassMethod,'DriftPass');
    wigglers=cellfun(iswiggler, ring);
    I2w=sum(cellfun(@wiggler_i2,ring(wigglers)));
    % Additional radiation
    I2x=sum(atgetfieldvalues(ring(atgetcells(ring,'I2')),'I2'));
    I2=nbper*(I2d+I2w+I2x);                                     % [m-1]
    U0=Cgamma/2/pi*(energy*1.e-9)^4*I2*1e9;                     % [eV]
end

    function i2=wiggler_i2(elem)
        rhoinv=elem.Bmax/Brho;
        coefh=elem.By(2,:);
        coefv=elem.Bx(2,:);
        i2=elem.Length*(coefh*coefh'+coefv*coefv')*rhoinv*rhoinv/2;
    end

end
