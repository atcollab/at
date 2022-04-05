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
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dp',dp)
%   Set the cavity frequency to the nominal value for the specified dp
%
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dct',dct)
%   Set the cavity frequency to the nominal value for the specified dct
%
%NEWRING=ATSETCAVITY(RING,...,'Voltage',VOLTAGE,...)
%   Set the total voltage (all cells) [V]
%   The voltage of each cavity is VOLTAGE / N_CAVITIES / PERIODICITY
%
%NEWRING=ATSETCAVITY(RING,...,'HarmNumber',h,...)
%   Set the harmonic number
%
%NEWRING=ATSETCAVITY(RING,...,'TimeLag',TIMELAG,...)
%   Set the time lag [m]
%
%NEWRING=ATSETCAVITY(RING,...,'cavpts',CAVPTS)
%   CAVPTS is the location of RF cavities. This allows to ignore harmonic
%   cavities. The default is to use all cavities
%
%  NOTES
%  1. In this mode, the radiation state of the lattice is not modified.
%  2. When dp is specified, the RF frequency is computed with the
%     slip factor, so that the resulting dp may slightly differ from the
%     specified value.
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

[props,idx]=atGetRingProperties(ring);
[cavpts,varargs]=getcavptsarg(varargin,ring,props);
[frequency,varargs]=getoption(varargs, 'Frequency', []);
[harmnumber,varargs]=getoption(varargs, 'HarmNumber',[]);
[vring,varargs]=getoption(varargs, 'Voltage', []);
[timelag,varargs]=getoption(varargs, 'TimeLag', []);
[dp,varargs]=getoption(varargs,'dp',NaN);
[dct,varargs]=getoption(varargs,'dct',NaN);

ncells=props.Periodicity;
cavities=ring(cavpts);
ncavs=length(cavities);
if ncavs == 0
    error('AT:NoCavity', 'No cavity found in the lattice');
end

if isempty(varargs)             % New syntax
    if ~isempty(frequency)
        lcell=findspos(ring,length(ring)+1);
        gamma0=props.Energy/props.Particle.rest_energy;
        beta0=sqrt(1-1/gamma0/gamma0);
        frev=beta0*CLIGHT/lcell/ncells;
        if (ischar(frequency) || isstring(frequency)) && strcmp(frequency, 'nominal')
            if isfinite(dct)
                frev=frev * lcell/(lcell+dct);
            elseif isfinite(dp)
                [~,ringrad]=check_radiation(ring,false,'force');
                etac=1/gamma0^2 - mcf(ringrad);
                frev=frev + frev*etac*dp;
            end
            frequency = frev * props_harmnumber(harmnumber,props);
        else
            harmnumber=round(frequency/frev);
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
    [vcell,radflag,harmnumber]=deal(varargin{:});
    harmcell=harmnumber/ncells;
    lcell=findspos(ring,length(ring)+1);
    gamma0=props.Energy/props.Particle.rest_energy;
    beta0=sqrt(1-1/gamma0/gamma0);
    frequency = (beta0*CLIGHT/lcell)*harmcell;

    %now set cavity frequencies, Harmonic Number and RF Voltage
    cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    cavities=atsetfieldvalues(cavities, 'HarmNumber', harmnumber);
    cavities=atsetfieldvalues(cavities, 'Voltage', vcell/ncavs);
    ring(cavpts)=cavities;

    %now set phaselags in cavities
    if radflag
        U0=atgetU0(ring);
        timelag= (lcell/(2*pi*harmcell))*asin(U0/vcell/ncells);
        ring=atradon(ring);  % set radiation on. nothing if radiation is already on
    else
        ring=atradoff(ring,'RFCavityPass');  % set radiation off. nothing if radiation is already off
        timelag=0;
    end
    ring=atsetfieldvalues(ring, cavpts, 'TimeLag', timelag);
end

% Update the ring properties if HarmNumber has changed
if ~(isempty(idx) || isempty(harmnumber) || ...
    (isfield(props, 'HarmNumber') && harmnumber == props.HarmNumber))
    ring=atSetRingProperties(ring,'HarmNumber',harmnumber);
end

    function h=props_harmnumber(h, props)
        if isempty(h)
            if isfield(props,'HarmNumber') && props.HarmNumber ~= 0
                h=props.HarmNumber;
            else
                error('AT:NoCavity','Harmonic number is unknown')
            end
        end
    end
end
