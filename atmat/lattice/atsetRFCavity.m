function newring = atsetRFCavity(ring, rfv, radflag, HarmNumber, DeltaFreq )
%  ATSETRFCAVITY - sets the RF Cavity with the passmethod RFCavityPass.
%  RFCavityPass allows to change the energy of the beam changing the
%  frequency of the cavity.
%
%   newring = atSetRFCavity(ring, rfv, radflag, HarmNumber, DeltaFreq)
%   sets the synchronous phase of the cavity, the voltage, the harmonic
%   number and the PassMethod. 
%
%  INPUTS
%  1. ring       - Ring structure
%  2. rfv        - RF-voltage in volts
%  3. radflag    - 0/1: activat/desactivate radiation (atradon/atradoff)
%  4. HarmNumber - Harmonic number
%  5. DeltaFreq  - Frequency shift in Hz
%
%  OUTPUTS
%  1. newring - update ring structure with nw RF parameters
%
%  NOTES
%  1. All the N cavities will have a voltage rfv/N
%  2. radflag says whether or not we want radiation on, which affects
%     synchronous phase. If radflag is 0, the function calls atradoff, 
%     if it is 1, it calls atradon.
%  3. Cavities in the ring must have the Class RFCavity.
%  4. Normally DeltaFreq should be 0, it's different from 0 when you want to
%     simulate a change of energy changing the RF frequency. DeltaFreq is in
%      Hz.
%  5. Does not work well for misaligned cavity 
%
%  EXAMPLES
%
%   1. normal use:
%   newring = atsetRFCavity( ring, 6e6, 1, 992, 0 )
%
%   2. for off-energy simulations:
%   newring = atsetRFCavity( ring, 6e6, 1, 992, 100 )
%
%  See also atsetcavity, atradon, atradoff, atgetU0

% Speed of light
CLIGHT = PhysConstant.speed_of_light_in_vacuum.value ;

newring = ring;

% Indices of all cavity
indrfc = findcells(ring,'Class','RFCavity');
beta0  = 1;
% get enegery losss per turn
U0     = atgetU0(ring); %% not ok if misaligned elem!

% Compute circumference
circ = findspos(ring,length(ring)+1);

% Set cavity pass
newring = setcellstruct(newring,'PassMethod',indrfc,'RFCavityPass');

% Set main frequency
freq0 = (CLIGHT*beta0/circ)*HarmNumber;
% Add offset
freq = freq0+DeltaFreq;

%now set cavity frequencies, Harmonic Number and RF Voltage
for j=indrfc
    newring{j}.Frequency=freq;
    newring{j}.HarmNumber=HarmNumber;
    newring{j}.Voltage=rfv/length(indrfc);
end

%now set phaselags in cavities
if radflag
    timelag= (circ/(2*pi*HarmNumber))*asin(U0/(rfv));
    newring=atradon(newring,'RFCavityPass','auto','auto');  % set radiation on. nothing if radiation is already on
    newring=setcellstruct(newring,'PassMethod',indrfc,'RFCavityPass');
else
    newring=atradoff(newring,'RFCavityPass','auto','auto');  % set radiation off. nothing if radiation is already off
    newring=setcellstruct(newring,'PassMethod',indrfc,'RFCavityPass');
    timelag=0;
end

newring = setcellstruct(newring,'TimeLag',indrfc,timelag);
end
