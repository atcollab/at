function ring = atsetcavity(ring,varargin)
%ATSECAVITY Set the cavity parameters
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
%NEWRING=ATSETCAVITY(RING,...,'HarmNumber',HARM_NUMBER,...)
%   Set the cavity harmonic number of the full ring (all cells)
%   The harmonic number of each cavity is HARM_NUMBER / N_CELLS
%
%NEWRING=ATSETCAVITY(RING,...,'Voltage',VOLTAGE,...)
%   Set the total voltage (all cells) [V]
%   The voltage of each cavity is VOLTAGE / N_CAVITIES / N_CELLS
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
%  RFV          RF voltage [V]
%  RADFLAG      0/1: activate/desactivate radiation (atradon/atradoff)
%  HARMNUMBER 	Harmonic number
%
%  NOTES
%  1. This mode is deprecated and should be replaced by
%       RING=ATSETCAVITY(RING,'Frequency','nominal',...
%           'HarmNumber',HARM_NUMBER, 'Voltage',RFV)
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

[refpts,varargs]=getoption(varargin, 'refpts', atgetcells(ring, 'Frequency'));
[frequency,varargs]=getoption(varargs, 'Frequency', []);
[voltage,varargs]=getoption(varargs, 'Voltage', []);
[harmnumber,varargs]=getoption(varargs, 'HarmNumber', []);
[timelag,varargs]=getoption(varargs, 'TimeLag', []);

if islogical(refpts)
    refpts=find(refpts);
end

[~,ncells]=atenergy(ring);
cavities=ring(refpts);
ncavs=length(cavities);

if isempty(varargs)             % New syntax
    if ncavs == 0
        error('AT:NoCavity', 'No cavity found in the lattice');
    end
    if ~isempty(harmnumber)
        cavities=atsetfieldvalues(cavities, 'HarmNumber', harmnumber/ncells);
    end
    if ~isempty(frequency)
        if (ischar(frequency) || isstring(frequency)) && strcmp(frequency, 'nominal')
            if isempty(harmnumber)
                harmnumber=ncells*getfield(cavities,'HarmNumber');
            end
            circ=ncells*findspos(ring,length(ring)+1);
            frequency = (CLIGHT/circ)*harmnumber;
        end
        cavities=atsetfieldvalues(cavities, 'Frequency', frequency);
    end
    if ~isempty(voltage)
        cavities=atsetfieldvalues(cavities, 'Voltage', voltage/ncells/ncavs);
    end
    if ~isempty(timelag)
        cavities=atsetfieldvalues(cavities, 'TimeLag', timelag);
    end
    ring(refpts)=cavities;
else                            % Old syntax, for compatibility
    % me_EV=510998.928;
    
    % gamma0=E0/me_EV;
    % beta0=sqrt(gamma0^2-1)/gamma0;
    [rfv,radflag,HarmNumber]=deal(varargin{:});
    L=findspos(ring,length(ring)+1);
    circ=L*ncells;
    %freq=(beta0*clight/circ)*HarmNumber;
    freq=(CLIGHT/circ)*HarmNumber;
    
    %now set cavity frequencies, Harmonic Number and RF Voltage
    cavities=atsetfieldvalues(cavities, 'Frequency', freq);
    cavities=atsetfieldvalues(cavities, 'HarmNumber', HarmNumber);
    cavities=atsetfieldvalues(cavities, 'Voltage', rfv/ncavs);
    ring(refpts)=cavities;
    
    %now set phaselags in cavities
    if radflag
        warning('AT:CavityTimeLag',...
            ['\nThis function modifies the time reference\n',...
            'This should be avoided, you have been warned\n']);
        U0=atgetU0(ring);
        timelag= (circ/(2*pi*HarmNumber))*asin(U0/(rfv));
        ring=atradon(ring);  % set radiation on. nothing if radiation is already on
    else
        ring=atradoff(ring,'CavityPass');  % set radiation off. nothing if radiation is already off
        timelag=0;
    end
    ring=atsetfieldvalues(ring, refpts, 'TimeLag', timelag);
end

    function values=getfield(ring,attr)
        values = unique(atgetfieldvalues(ring,attr));
        if length(values) > 1
            error('AT:IncompatibleValues','%s not equal for all cavities', attr);
        end
    end
        
end
