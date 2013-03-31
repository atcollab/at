% fit a bump using correctors 
clear all

load('dba.mat','RING');
addpath(fullfile(pwd,'..'))

%correctors and BPM
C=atcorrector('C',0,0);
M=atmarker('BPM');

% get one cell and add elements
arc=[{M};RING(1:18);RING(128:end)];

indq=findcells(arc,'Class','Quadrupole');
for iq=2:2:length(indq)
    arc=[arc(1:indq(iq)-1);M;C;arc(indq(iq):end)];
    indq=findcells(arc,'Class','Quadrupole');
end

% build variables
hcor=findcells(arc,'FamName','C');

Variab=atVariableBuilder(arc,...
    {[hcor(1), hcor(end)],[hcor(2),hcor(end-1)]},...
    {{'KickAngle'}});

% build constraints
bpm=findcells(arc,'FamName','BPM');       
          
c1=atlinconstraint(...
    [bpm(1)],...
    {{'ClosedOrbit',{1}},{'ClosedOrbit',{2}}},...
    [1e-3,0],...
    [1e-3,0],...
    [1e-2 1e-2]);

c2=atlinconstraint(...
    [bpm(2:end-1)],{{'ClosedOrbit',{1}}},0,0,1e-2); %#ok<*NBRAK>

c=[c1,c2];

% perform matching
arc_bump=atmatch(arc,Variab,c,10^-15,1000,3,@lsqnonlin);%'fminsearch',3);%
figure;atplot(arc_bump,@plClosedOrbit);


