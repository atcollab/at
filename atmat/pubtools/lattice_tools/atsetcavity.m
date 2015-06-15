function newring = atsetcavity(ring,rfv, radflag,HarmNumber)
%newring = atsetcavity(ring,rfv, radflag,HarmNumber)
%sets the synchronous phase of the cavity assuming radiation is turned on
%radflag says whether or not we want radiation on, which affects
%synchronous phase.
%also sets the rf voltage and Harmonic number
%also sets the rf frequency.

clight=2.99792458e8 ;
me_EV=510998.928;

newring = ring;

indrfc=findcells(ring,'Class','RFCavity');

E0=ring{indrfc(1)}.Energy;
gamma0=E0/me_EV;
beta0=sqrt(gamma0^2-1)/gamma0;

U0=atgetU0(ring);

%find circumference
nc=2*pi/sum(getcellstruct(ring,'BendingAngle',findcells(ring,'BendingAngle')));
ncells=round(nc);
if ~isfinite(nc)
    warning('AT:WrongNumberOfCells','No bending in the cell, ncells set to 1');
    ncells=1;
elseif abs(nc-ncells) > 1.e-4
    warning('AT:WrongNumberOfCells','non integer number of cells: ncells = %g',nc);
end

L=findspos(ring,length(ring)+1);
circ=L*ncells;
%freq=(beta0*clight/circ)*HarmNumber;
freq=(clight/circ)*HarmNumber;

%now set cavity frequencies, Harmonic Number and RF Voltage
for j=indrfc
    newring{j}.Frequency=freq;
    newring{j}.HarmNumber=HarmNumber;
    newring{j}.Voltage=rfv/length(indrfc);
end

%now set phaselags in cavities
if radflag
    timelag= (circ/(2*pi*HarmNumber))*asin(U0/(rfv));
    newring=atradon(newring);  % set radiation on. nothing if radiation is already on
    
else
    newring=atradoff(newring);  % set radiation off. nothing if radiation is already off
    timelag=0;
end

newring=setcellstruct(newring,'TimeLag',indrfc,timelag);
end
