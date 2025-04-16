function ring = atsetcavity(ring,varargin)
%ATSECAVITY Set the parameters of RF cavities
%
%ATSETCAVITY may be used in two modes:
%
%Upgrade mode
%===================================================
% By default, ATSETCAVITY will act on the "main" cavities: they are defined by the
% cavpts ring property, or if absent by cavities at the lowest frequency.
%
%NEWRING=ATSETCAVITY(RING,...,'Frequency',FREQUENCY,...)
%   Set the cavity frequency [Hz]. FREQUENCY is a scalar or an array as
%   long as the list of selected cavities
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
%   Set the total voltage (all cells) [V]. VOLTAGE will be distributed over the
%   cells: CELL_VOLTAGE = VOLTAGE / PERIODICITY.
%   Then if CELL_VOLTAGE is a scalar, it will be equally shared among the
%   selected cavities. Otherwise it is an array as long as the list of
%   selected cavities.
%
%NEWRING=ATSETCAVITY(RING,...,'HarmNumber',H,...)
%   Set the harmonic number. H is a scalar or an array as
%   long as the list of selected cavities
%
%NEWRING=ATSETCAVITY(RING,...,'TimeLag',TIMELAG,...)
%   Set the time lag [m], . TIMELAG is a scalar or an array as
%   long as the list of selected cavities
%
%NEWRING=ATSETCAVITY(RING,...,'cavpts',CAVPTS)
%   CAVPTS is the location of the selected RF cavities. The default is to act on the
%   "main" cavities: they are defined by the cavpts ring property, or if absent by
%   cavities at the lowest frequency.
%
%  NOTES
%  1. In this mode, the radiation state of the lattice is not modified.
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
%       RING=ATSETCAVITY(RING,'Frequency','nominal','HarmNumber',HARM_NUMBER, 'Voltage',RFV)
%       RING=atSetCavityPhase(RING) (optional)
%       RING=atenable_6d(RING)      (optional)
%  2. All the N cavities will have a voltage RFV/N
%  3. sets the synchronous phase of the cavity assuming radiation is turned
%     on radflag says whether or not we want radiation on, which affects
%     synchronous phase.
%
%  See also atSetCavityPhase, atsetRFcavity, atenable_6d, atdisable_6d, atgetU0

% Speed of light
CLIGHT=PhysConstant.speed_of_light_in_vacuum.value;

[ncells,cell_h,beta0,maincavs]=atGetRingProperties(ring,'Periodicity','cell_harmnumber','beta','cavpts');

[cavpts,varargs]=getoption(varargin,'cavpts',maincavs);
[frequency,varargs]=getoption(varargs, 'Frequency', []);
[harmnumber,varargs]=getoption(varargs, 'HarmNumber',[]);
[vring,varargs]=getoption(varargs, 'Voltage', []);
[vcell,varargs]=getoption(varargs, 'cell_voltage', []);
[timelag,varargs]=getoption(varargs, 'TimeLag', []);
[dp,varargs]=getoption(varargs,'dp',NaN);
[dct,varargs]=getoption(varargs,'dct',NaN);
[df,varargs]=getoption(varargs,'df',NaN);

cavities=ring(cavpts);
ncavs=length(cavities);
if ncavs == 0
    error('AT:NoCavity', 'No cavity found in the lattice');
end

if isempty(varargs)             % New syntax
    if ~isempty(harmnumber)
        cavities=atsetfieldvalues(cavities, 'HarmNumber', harmnumber/ncells);
    end
    if ~isempty(frequency)
        if ischar(frequency) || isstring(frequency)
            if strcmp(frequency,'nominal')
                lcell=findspos(ring,length(ring)+1);
                frev=beta0*CLIGHT/lcell;
                hh = cellfun(@getcavh, cavities);
                if isfinite(df)
                    frev = frev + df/min(hh);
                elseif isfinite(dct)
                    frev=frev * lcell/(lcell+dct);
                elseif isfinite(dp)
                    % Find the path lengthening for dp
                    [~,rnorad]=check_radiation(ring,false,'force');
                    [~,orbitin]=findorbit4(rnorad,dp, 'strict', -1);
                    orbitout=ringpass(rnorad,orbitin);
                    dct=orbitout(6);
                    frev=frev * lcell/(lcell+dct);
                end
                frequency = hh * frev;
            else
                error('AT:Frequency', 'Wrong frequency specifiation: ''%s''',frequency);
            end
        end
        cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    end
    if ~isempty(vring)
        if numel(vring) ~= ncavs, vring=vring/ncavs; end
        cavities=atsetfieldvalues(cavities, 'Voltage', vring/ncells);
    end
    if ~isempty(vcell)
        if numel(vcell) ~= ncavs, vcell=vcell/ncavs; end
        cavities=atsetfieldvalues(cavities, 'Voltage', vcell);
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
if ~isempty(idx)
    h=unique(atgetfieldvalues(ring,maincavs,'HarmNumber'));
    if h(1) ~= cell_h 
        ring=atSetRingProperties(ring,'cell_harmnumber',h(1));
    end
end

    function h=getcavh(cav,frev)
        if isfield(cav,'HarmNumber')
            h=cav.HarmNumber;
        else
            h=round(cav.Frequency/frev);
        end
    end
end
