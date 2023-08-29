function ring = atsetcavity(ring,varargin)
%ATSECAVITY Set the parameters of RF cavities
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
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dp',DP)
%   Set the cavity frequency to the nominal value for the specified dp
%
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','dct',DCT)
%   Set the cavity frequency to the nominal value for the specified dct
%
%NEWRING=ATSETCAVITY(RING,...,'Frequency','nominal','df',DF)
%   Set the cavity frequency to the nominal value + df
%
%NEWRING=ATSETCAVITY(RING,...,'Voltage',VOLTAGE,...)
%   Set the total voltage (all cells) [V]
%   The voltage of each cavity is VOLTAGE / N_CAVITIES / PERIODICITY
%
%NEWRING=ATSETCAVITY(RING,...,'HarmNumber',H,...)
%   Set the harmonic number
%
%NEWRING=ATSETCAVITY(RING,...,'TimeLag',TIMELAG,...)
%   Set the time lag [m]
%
%NEWRING=ATSETCAVITY(RING,...,'cavpts',CAVPTS)
%   CAVPTS is the location of RF cavities. The default is to act on the
%   "main" cavities: if there is is a 'cavpts' ring property. it defined
%   the main cavities, otherwise the main cavities are cavities with the
%   lowest frequency
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
%  RFV          RF voltage (full ring) [V]
%  RADFLAG      0/1: activate/desactivate radiation (atradon/atradoff)
%  HARMNUMBER 	Harmonic number (full ring)
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
%  See also atSetCavityPhase, atsetRFcavity, atenable_6d, atdisable_6d, atgetU0

% Speed of light
CLIGHT=PhysConstant.speed_of_light_in_vacuum.value;

[cavpts,varargs]=getoption(varargin,'cavpts',[]);
[frequency,varargs]=getoption(varargs, 'Frequency', []);
[harmnumber,varargs]=getoption(varargs, 'HarmNumber',[]);
[vring,varargs]=getoption(varargs, 'Voltage', []);
[vcell,varargs]=getoption(varargs, 'cell_voltage', []);
[timelag,varargs]=getoption(varargs, 'TimeLag', []);
[dp,varargs]=getoption(varargs,'dp',NaN);
[dct,varargs]=getoption(varargs,'dct',NaN);
[df,varargs]=getoption(varargs,'df',NaN);

if isempty(cavpts)
    [ncells,cell_h,beta0,cavpts]=atGetRingProperties(ring,'Periodicity','cell_harmnumber','beta','cavpts');
else
    [ncells,cell_h,beta0]=atGetRingProperties(ring,'Periodicity','cell_harmnumber','beta');
end
cavities=ring(cavpts);
ncavs=length(cavities);
if ncavs == 0
    error('AT:NoCavity', 'No cavity found in the lattice');
end

if isempty(varargs)             % New syntax
    if isempty(harmnumber)
        harmcell=[];
    else
        harmcell=harmnumber/ncells;
    end
    if ~isempty(frequency)
        lcell=findspos(ring,length(ring)+1);
        frev=beta0*CLIGHT/lcell;
        if (ischar(frequency) || isstring(frequency)) && strcmp(frequency, 'nominal')
            hh=props_harmnumber(harmcell,cell_h);
            if isfinite(df)
                frev = frev + df/hh;
            elseif isfinite(dct)
                frev=frev * lcell/(lcell+dct);
            elseif isfinite(dp)
                % Find the path lengthening for dp
                [~,rnorad]=check_radiation(ring,false,'force');
                [~,orbitin]=findorbit4(rnorad,dp);
                orbitout=ringpass(rnorad,orbitin);
                dct=orbitout(6);
                frev=frev * lcell/(lcell+dct);
            end
            frequency = hh * frev;
        else
            harmcell=round(frequency/frev);
        end
        cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    end
    if ~isempty(vring)
        cavities=atsetfieldvalues(cavities, 'Voltage', vring/ncells/ncavs);
    end
    if ~isempty(vcell)
        cavities=atsetfieldvalues(cavities, 'Voltage', vring/ncavs);
    end
    if ~isempty(timelag)
        cavities=atsetfieldvalues(cavities, 'TimeLag', timelag);
    end
    ring(cavpts)=cavities;
else                            % Old syntax, for compatibility
    [vring,radflag,harmnumber]=deal(varargin{:});
    harmcell=harmnumber/ncells;
    lcell=findspos(ring,length(ring)+1);
    frequency = (beta0*CLIGHT/lcell)*harmcell;

    %now set cavity frequencies, Harmonic Number and RF Voltage
    cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    cavities=atsetfieldvalues(cavities, 'HarmNumber', harmnumber);
    cavities=atsetfieldvalues(cavities, 'Voltage', vring/ncavs);
    ring(cavpts)=cavities;

    %now set phaselags in cavities
    if radflag
        U0=atgetU0(ring);
        timelag= (lcell/(2*pi*harmcell))*asin(U0/vring);
        ring=atenable_6d(ring);  % set radiation on. nothing if radiation is already on
    else
        ring=atdisable_6d(ring,'cavipass','RFCavityPass');  % set radiation off and turn on cavities
        timelag=0;
    end
    ring=atsetfieldvalues(ring, cavpts, 'TimeLag', timelag);
end

% Update the ring properties if HarmNumber has changed
idx=atlocateparam(ring);
if ~(isempty(idx) || isempty(harmcell) || ...
    (harmcell == cell_h))
    ring=atSetRingProperties(ring,'cell_harmnumber',harmcell);
end

    function cellharm=props_harmnumber(cellharm, cell_h)
        if isempty(cellharm)
            if ~isfinite(cell_h)
                error('AT:NoCavity','Harmonic number is unknown')
            end
            cellharm=cell_h;
        end
    end
end
