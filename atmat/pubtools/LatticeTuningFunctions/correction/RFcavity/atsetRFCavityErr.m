function [newring, inCOD]= atsetRFCavityErr( ring, rfv, radflag, HarmNumber, inCOD, DeltaFreq )
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
%   modified : defaulted DeltaFreq
%              removed U0 set timelag to 6th coordinate of findorbit6
%
%
%   see also: atsetcavity

clight=299792458 ;

newring = ring;
indrfc=findcells(ring,'Class','RFCavity');

beta0=1;

% circumference + trajectory length modification introduced by orbit in
% dipoles

% DLL=LatticeLengtheningCODDipole(ring,inCOD);
DLL=0;
circumference=findspos(ring,length(ring)+1)*(1+DLL);

newring=setcellstruct(newring,'PassMethod',indrfc,'RFCavityPass');

freq0=(clight*beta0/circumference)*HarmNumber;

%now set intial guess for cavity frequencies,
% set Harmonic Number and RF Voltage
for j=indrfc
    newring{j}.Frequency=freq0;
    newring{j}.HarmNumber=HarmNumber;
    newring{j}.Voltage=rfv/length(indrfc);
end

% set radiation pass methods ('auto' sets radiation on also in quadrupoles)
if radflag
    newring=atradon(newring,'RFCavityPass','auto','auto');
else
    newring=atradoff(newring,'RFCavityPass','auto','auto');
end

% default DeltaFreq to cancel energy deviation after one turn
if nargin<6
    disp('DeltaFreq not provided. Setting DeltaFreq to cancel energy variation after one turn');
    
    for i=1:5
        freq=atgetfieldvalues(newring,indrfc,'Frequency');
        tlag=atgetfieldvalues(newring,indrfc,'TimeLag');
        
        %alpha=mcfErr(newring,indrfc,inCOD,DLL);
        alpha=mcf(newring);
        orb = findorbit6(newring,indrfc,inCOD);
        if ~isnan(orb(1))
            df = alpha*orb(5,:)'.*freq;
        else
            df=zeros(size(freq));
            warning('findorbit6 failed.')
        end
        
        newring=atsetfieldvalues(newring,indrfc,'Frequency',freq+df);
        
        orb = findorbit6(newring,indrfc,inCOD);
        newring=atsetfieldvalues(newring,indrfc,'TimeLag',tlag-orb(6,:)');
    
        
    end
    
    freqct=atgetfieldvalues(newring(indrfc),'Frequency');
    DeltaFreq=freqct(1)-freq0;
    
    
    disp(['DeltaFreq: ' num2str(DeltaFreq) ' Hz']);
    
end

freq=freq0+DeltaFreq;

% set frequency
newring(indrfc)=atsetfieldvalues(newring(indrfc),'Frequency',freq);

% set timelag
orb = findorbit6(newring,indrfc,inCOD);
tlag=atgetfieldvalues(newring,indrfc,'TimeLag');
ntlag=tlag-orb(6,:)';

newring(indrfc)=atsetfieldvalues(newring(indrfc),'TimeLag',ntlag);

inCOD = findorbit6(newring,1,inCOD);

if nargin<7
    verbose=false;
end

if verbose
    disp(['DLL: ' num2str(circumference*DLL) ' m']);
    disp(['Delta RF frequency: ' num2str(DeltaFreq) ' Hz']);
    disp(['Time Lag: ' num2str(ntlag) ' m']);
end

end
