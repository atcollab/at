function newring = atsetRFCavity( ring, rfv, radflag, HarmNumber, DeltaFreq )
%  ATSETRFCAVITY sets the RF Cavity with the passmethod RFCavityPass.
%  RFCavityPass allows to change the energy of the beam changing the
%  frequency of the cavity.
%
%   newring = atSetRFCavity(ring, rfv, radflag, HarmNumber, DeltaFreq)
%   sets the synchronous phase of the cavity, the voltage, the harmonic
%   number and the PassMethod. All the N cavities will have a voltage rfv/N
%   radflag says whether or not we want radiation on, which affects
%   synchronous phase. If radflag is 0, the function calls atradoff, if it
%   is 1, it calls atradon.
%   Cavities in the ring must have the Class RFCavity.
%   Normally DeltaFreq should be 0, it's different from 0 when you want to
%   simulate a change of energy changing the RF frequency. DeltaFreq is in
%   Hz.
%
%   normal use:
%   newring = atsetRFCavity( ring, 6e6, 1, 992, 0 )
%   for off-energy simulations:
%   newring = atsetRFCavity( ring, 6e6, 1, 992, 100 )
%
%   see also: atsetcavity

clight=299792458 ;

newring = ring;
indrfc=findcells(ring,'Class','RFCavity');
beta0=1;
U0=atgetU0(ring); %% not ok if misaligned elem!

circ=findspos(ring,length(ring)+1);

newring=setcellstruct(newring,'PassMethod',indrfc,'RFCavityPass');

freq0=(clight*beta0/circ)*HarmNumber;

freq=freq0+DeltaFreq;

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
newring=setcellstruct(newring,'TimeLag',indrfc,timelag);
end
