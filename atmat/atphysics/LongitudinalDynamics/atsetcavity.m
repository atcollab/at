function ring = atsetcavity(ring,varargin)
%ATSECAVITY Set the cavity parameters
%
%WARNING: This function modifies the time reference,
%this should be avoided
%
%ATSETCAVITY may be used in two modes:
%
%Upgrade mode
%===================================================
%NEWRING=ATSETCAVITY(RING,...,'Frequency',FREQUENCY,...)
%   Set the cavity frequency [Hz].
%
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal',...)
%   Set the cavity frequency to the nominal value according to
%   circumference and harmonic number
%
%NEWRING=ATSETCAVITY(RING,...,'Voltage',VOLTAGE,...)
%   Set the total voltage (all cells) [V]
%   The voltage of each cavity is VOLTAGE / N_CAVITIES / PERIODICITY
%
%NEWRING=ATSETCAVITY(RING,...,'TimeLag',TIMELAG,...)
%   Set the time lag [m]
%
%NEWRING=ATSETCAVITY(RING,...,'refpts',CAVPTS)
%   CAVPTS is the location of RF cavities. This allows to ignore harmonic
%   cavities. The default is to use all cavities
%
%  NOTES
%  1. In this mode, the radiation state of the lattice is not modified.
%
%
%Compatibility mode
%===================================================
%NEWRING = ATSETCAVITY(RING,RFV,RADFLAG,HARM_NUMBER)
%  RING         Ring structure
%  RFV          RF voltage (per cell) [V]
%  RADFLAG      0/1: activate/desactivate radiation (atradon/atradoff)
%  HARMNUMBER 	Harmonic number (for 1 cell)
%
%  NOTES
%  1. This mode is deprecated and should be replaced by
%       RING=ATSETCAVITY(RING,'Frequency','nominal',...
%           'HarmNumber',HARM_NUMBER*PERIODICITY, 'Voltage',RFV/PERIODICITY)
%       RING=atSetCavityPhase(RING) (optional)
%       RING=atradon(RING)          (optional)
%  2. All the N cavities will have a voltage RFV/N
%  3. sets the synchronous phase of the cavity assuming radiation is turned
%     on radflag says whether or not we want radiation on, which affects
%     synchronous phase.
%
%  See also atSetCavityPhase, atsetRFcavity, atradon, atradoff, atgetU0

% Speed of light
CLIGHT=PhysConstant.speed_of_light_in_vacuum.value;
E_MASS=1.0E6*PhysConstant.electron_mass_energy_equivalent_in_MeV.value;

props=atGetRingProperties(ring);
try         % Look for cavities in the lattice properties
    cavpts=props.cavpts;
catch       % Take all cavities
    cavpts=atgetcells(ring,'Frequency');
end
[cavpts,varargs]=getoption(varargin, 'refpts', cavpts);
[frequency,varargs]=getoption(varargs, 'Frequency', []);
[vring,varargs]=getoption(varargs, 'Voltage', []);
[timelag,varargs]=getoption(varargs, 'TimeLag', []);
[dp,varargs]=getoption(varargs,'dp',NaN);
[dct,varargs]=getoption(varargs,'dct',NaN);

if islogical(cavpts)
    cavpts=find(cavpts);
end

ncells=props.Periodicity;
cavities=ring(cavpts);
ncavs=length(cavities);

if isempty(varargs)             % New syntax
    if ncavs == 0
        error('AT:NoCavity', 'No cavity found in the lattice');
    end
    if ~isempty(frequency)
        if (ischar(frequency) || isstring(frequency)) && strcmp(frequency, 'nominal')
            gamma0=props.Energy/E_MASS;
            beta0=sqrt(1-1/gamma0/gamma0);
            lcell=ncells*findspos(ring,length(ring)+1);
%           frev=beta0*CLIGHT/lcell;
            frev=CLIGHT/lcell;
            if isfinite(dct)
                frev=frev - frev*frev/CLIGHT*ncells*dct;
            elseif isfinite(dp)
                [~,ringrad]=check_radiation(ring,false,'force');
                etac=1/gamma0^2 - mcf(ringrad);
                frev=frev + frev*etac*dp;
            end
            frequency = frev*props.HarmNumber;
        end
        cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    end
    if ~isempty(vring)
        cavities=atsetfieldvalues(cavities, 'Voltage', vring/ncells/ncavs);
    end
    if ~isempty(timelag)
        cavities=atsetfieldvalues(cavities, 'TimeLag', timelag);
    end
    ring(cavpts)=cavities;
else                            % Old syntax, for compatibility
    % gamma0=props.Energy/E_MASS;
    % beta0=sqrt(gamma0^2-1)/gamma0;
    [vcell,radflag,harmcell]=deal(varargin{:});
    lcell=findspos(ring,length(ring)+1);
%   frequency = (beta0*CLIGHT/circ)*harmnumber;
    frequency = (CLIGHT/lcell)*harmcell;
    
    %now set cavity frequencies, Harmonic Number and RF Voltage
    cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    cavities=atsetfieldvalues(cavities, 'HarmNumber', harmcell);
    cavities=atsetfieldvalues(cavities, 'Voltage', vcell/ncavs);
    ring(cavpts)=cavities;
    
    %now set phaselags in cavities
    if radflag
        U0=atgetU0(ring);
        timelag= (lcell/(2*pi*harmcell))*asin(U0/vcell/ncells);
        ring=atradon(ring);  % set radiation on. nothing if radiation is already on
    else
        ring=atradoff(ring,'CavityPass');  % set radiation off. nothing if radiation is already off
        timelag=0;
    end
    ring=atsetfieldvalues(ring, cavpts, 'TimeLag', timelag);
end

    function values=getfield(ring,attr)
        values = unique(atgetfieldvalues(ring,attr));
        if length(values) > 1
            error('AT:IncompatibleValues','%s not equal for all cavities', attr);
        end
    end
        
end
