function newring = atsetcavity(ring,rfv,radflag,HarmNumber)
%ATSECAVITY Sets the synchronous phase of the cavity assuming radiation
%
%  newring = atsetcavity(ring,rfv, radflag,HarmNumber)
%
%  INPUTS
%  1. ring       - Ring structure
%  2. rfv        - RF-voltage in volts
%  3. radflag    - 0/1: activat/desactivate radiation (atradon/atradoff)
%  4. HarmNumber - Harmonic number
%
%  OUTPUTS
%  1. newring - Updated ring structure with nw RF parameters
%
%  NOTES
%  1. All the N cavities will have a voltage rfv/N
%  2. sets the synchronous phase of the cavity assuming radiation is turned
%     on radflag says whether or not we want radiation on, which affects
%     synchronous phase.
%
%  See also atsetRFcavity, atradon, atradoff, atgetU0


% Speed of light
CLIGHT=PhysConstant.speed_of_light_in_vacuum.value;
% me_EV=510998.928;

newring = ring;

[~,ncells]=atenergy(ring);
% gamma0=E0/me_EV;
% beta0=sqrt(gamma0^2-1)/gamma0;

L=findspos(ring,length(ring)+1);
circ=L*ncells;
%freq=(beta0*clight/circ)*HarmNumber;
freq=(CLIGHT/circ)*HarmNumber;

%now set cavity frequencies, Harmonic Number and RF Voltage
indrfc=findcells(ring,'Class','RFCavity');
for j=indrfc
    newring{j}.Frequency=freq;
    newring{j}.HarmNumber=HarmNumber;
    newring{j}.Voltage=rfv/length(indrfc);
end

%now set phaselags in cavities
if radflag
    U0=atgetU0(ring);
    timelag= (circ/(2*pi*HarmNumber))*asin(U0/(rfv));
    newring=atradon(newring);  % set radiation on. nothing if radiation is already on    
else
    newring=atradoff(newring,'CavityPass');  % set radiation off. nothing if radiation is already off
    timelag=0;
end

newring=setcellstruct(newring,'TimeLag',indrfc,timelag);
end
