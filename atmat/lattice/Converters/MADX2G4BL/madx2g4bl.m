% function [outtext]=madx2g4bl(P_0,particle,folder)
% tansform madx file into G4BL input file.
%
% Simone Maria Liuzzo PhD@LNF 11-jul-2011

format longe

StepLength=10; % [mm] default is 10mm

Nevents=1; % e+ events and e- events.
disp(['writing G4BL file to simulate: ' num2str(Nevents) ' e+ and ' num2str(Nevents) ' e-'])
det=0;
SR=0;

TotLength=1258358.149; %[mm]

%beamstart=TotLength/2;
%beamstart=581.981242*1000-140; % first QDI center

beamstart=0;

% beams shapes

% diversi da zero muore!!! sestupoli??? KL? BRHO?
sX=1*1e-3; %[mm]
sY=1*1e-3; %[mm]
sXp=1*1e-6; %[rad]
sYp=1*1e-6; %[rad]

killval=1;

Xang=0.06600162900000006; % rad crossing angle


%%%% HER %%%%%
%%%%%%%%%%%%%%

folder='SBHERV12DATA';
P_0=6.699999981; % GeV/c
particle='e+';

sxt=1;
quadH=1;  % go around in the other direction.

Brho=3.3356409519815205*P_0; % positron and opposite rotation


bigmondo=[...
    ...'box mondo height=10000 width=3000000 length=3000000 material=Vacuum \n'...
    ];

posbeam=[ ' beam gaussian particle=' particle ...
    ' beamZ=' num2str(beamstart) ...
    ' beamX=' num2str(+1*1.532845188e-07,15) ...
    ' beamXp=' num2str(+1*1.683548689e-07,15) ...
    ' rotation=Y' num2str(0*180)... rotate to head back
    '  nEvents=' num2str(Nevents) ' '... % nEvents=10
    '	sigmaX=' num2str(sX,15) ...
    ' sigmaY=' num2str(sY,15)...
    ' sigmaXp=' num2str(sXp,15)...
    ' sigmaYp=' num2str(sYp,15) ' '...
    '	meanMomentum=' num2str(P_0*1000,15) ' sigmaP=0'...
    ' meanT=0 sigmaT=0 \n'...
    'trackcolor reference=1,1,1 \n'...
    ];

head=['physics QGSP synchrotronRadiation=' num2str(SR) '\n'...
    'start x=0 y=0 z=' num2str(beamstart) ' radiusCut=0.001 ring=1\n'...
    'reference particle=' particle ' referenceMomentum=' num2str(P_0*1000,15) ' noEloss=1'...
    ' beamXp=0 beamYp=0 '...
    ' tuneMomentum=' num2str(P_0,15) '\n'...
    'box BeamVis width=200.0 height=200.0 length=0.1 material=Vacuum color=1,0,0 kill=1 \n'...
    'trace nTrace=' num2str(2*Nevents) ' oneNTuple=1 format=ascii \n'...
    ' param eventTimeLimit=60*8\n'... %  8 minutes maximum time for one simulation
    ' trackcuts keep=e+,e-,gamma maxTime=1000000*40\n' ... %ns --> 4micros
    ];
% x damping time= 0.26700741E-01 => 15% =~ 40ms = 40*10^6 ns

% get pos K and L for Quad, sext and bend
Qfile=fopen([folder '/quad']);
Sfile=fopen([folder '/sext']);
Bfile=fopen([folder '/bend']);
% get pos K and L for Quad, sext and bend -- sbend and not rbend.
%Qfile=fopen('DIAMDATA/quad');
%Sfile=fopen('DIAMDATA/sext');
%Bfile=fopen('DIAMDATA/bend');

% name s K L
Q=textscan(Qfile,'%s %f %f %f','Headerlines',47);
B=textscan(Bfile,'%s %f %f %f','Headerlines',47);
S=textscan(Sfile,'%s %f %f %f','Headerlines',47);

defQ= ['virtualdetector Det radius=10.0 length=' num2str(2*StepLength)...
    ' color=0,.2,.0 referenceParticle=1 format=rootExtended\n'];
defS='';
defB='';
ord=1;
plac={};
defined={};
d=1;

for i=1:length(Q{1})
    name=(Q{1}(i));
    name{1}(1)='';
    name{1}(end)='';
    Q{1}(i)=name;
end
for i=1:length(S{1})
    name=S{1}(i);
    name{1}(1)='';
    name{1}(end)='';
    S{1}(i)=name;
end
for i=1:length(B{1})
    name=B{1}(i);
    name{1}(1)='';
    name{1}(end)='';
    B{1}(i)=name;
end

% small final focus quad
ffQuad{1}=char(Q{1}(1));%{'QD0AL'};
ffQuad{2}=char(Q{1}(2));%{'QD0BL'};
ffQuad{3}=char(Q{1}(3));%{'QD0AR'};
ffQuad{4}=char(Q{1}(4));%{'QD0AR'};
ffQuad{5}=char(Q{1}(end));%{'QD0BR'};
ffQuad{6}=char(Q{1}(end-1));%{'QDPMR'};
ffQuad{7}=char(Q{1}(end-2));%{'QDPML'};
ffQuad{8}=char(Q{1}(end-3));%{'QD0BR'};


plac{ord}=[' cornerarc' ' z=0+' num2str(0) '+0+' num2str(0,15) ...
    ' angle=' num2str(-Xang/2*180/pi,15)  ' centerRadius='  num2str(1,15) '\n'...
    ];
ord=ord+1;

for i=1:length(Q{1})
    %genericquad QDpi fieldLength=280 apertureRadius=101.5 ironRadius=381 \
    %	ironLength=250 ironColor=0.6,.0,.6 kill=' num2str(killval) ' gradient=$KQDpi
    %
    if sum(strcmp(char(Q{1}(i)),defined))==0;
        inrad=60;
        outrad=70;
        
        if sum(strcmp(char(Q{1}(i)),ffQuad))~=0
            inrad=12;
            outrad=15;
            tempSL=StepLength;
            StepLength=0.01;
            if strcmp(char(Q{1}(i)),ffQuad{1}) || strcmp(char(Q{1}(i)),ffQuad{5}) %% P.M.
                inrad=8;
                outrad=10;
                
            end
            
        end
        defQ=[ defQ ...
            ' genericquad ' Q{1}(i) ' '...
            ' fieldLength=' num2str(Q{4}(i)*1000,15) ... % mm
            ' apertureRadius=' num2str(inrad,15) ' '...
            ' ironRadius=' num2str(outrad,15) ' '...
            ' ironLength=' num2str(Q{4}(i)*1000,15) ' '...
            ' ironColor=0/255,191/255,255/255 kill=' num2str(killval) ' '...
            ' fringe=0'...
            ' gradient=' num2str(quadH*Q{3}(i)*Brho/Q{4}(i),15) ' '... K
            ' maxStep=' num2str(StepLength/2)  ...
            '\n']; %#ok<*AGROW>
        
        % add to already defined
        defined{d}=char(Q{1}(i));
        d=d+1;
    end
    
    StepLength=tempSL;
    
    %defQ=[defQ 'virtualdetector Det radius=400.0 color=0,.1,.1'];
    
    if (ord/2-floor(ord/2))==(det-1) % set 0 to place a detector in every quadrupole.
        plac{ord}=['place ' Q{1}(i) ' parent=  z=' num2str(Q{2}(i)*1000,15)...
            ' coordinates=centerline'...
            '\n'...
            'place ' 'Det' ...
            ' rename=DET# '...
            ' parent=  z=' num2str((Q{2}(i))*1000)...
            ' coordinates=centerline '...
            '\n']; %#ok<*SAGROW>
    else
        plac{ord}=['place ' Q{1}(i) ' parent=  z=' num2str(Q{2}(i)*1000,15)...
            ' coordinates=centerline'...
            '\n']; %#ok<*SAGROW>
    end
    %      if sum(strcmp(Q{1}(i),{'"QD0AL"';'"QD0BL"';'"QD0AR"';'"QD0BR"'}))>0
    %
    %      plac{ord}=['place ' Q{1}(i) ' parent=  z=' num2str(Q{2}(i)*1000,15)...
    %      ' coordinates=centerline'...
    %      ' rotation=Y' num2str(-0.066*180/pi,15) ' '...
    %      '\n'];
    %     end
    ord=ord+1;
end

for i=1:length(S{1})
    %multipole SF1 fieldLength=150 apertureRadius=101.5 ironRadius=381 \
    %	ironLength=140 ironColor=0,.6,.6 kill=' num2str(killval) ' sextupole=$KSF1
    
    if sum(strcmp(char(S{1}(i)),defined))==0;
        inrad=80;
        outrad=120;
        
        defS=[ defS ...
            ' multipole ' S{1}(i) ' '...
            ' fieldLength=' num2str(S{4}(i)*1000,15) ... % mm
            ' apertureRadius=' num2str(inrad,15) ' '...
            ' ironRadius=' num2str(outrad,15) ' '...
            ' ironLength=' num2str(S{4}(i)*1000,15) ' '...
            ' ironColor=173/255,255/255,47/255 kill=' num2str(killval) ' '...
            ' fringe=0'...
            ' maxStep=' num2str(StepLength/2)  ...
            ' sextupole=' num2str(sxt*S{3}(i)*Brho/S{4}(i),15) ' '... K2
            '\n'];
        % add to already defined
        defined{d}=char(S{1}(i));
        d=d+1;
    end
    
    
    plac{ord}=['place ' S{1}(i) ' parent=  z=' num2str(S{2}(i)*1000,15)...
        ' coordinates=centerline '...
        '\n']; %#ok<*SAGROW>
    ord=ord+1;
    
end

for i=1:length(B{1})
    
    %genericbend BHer fieldWidth=1660 fieldHeight=200 fieldLength=5400 \
    %	ironColor=1,0,0 ironWidth=1828 ironHeight=1320 ironLength=5000 \
    % kill=' num2str(killval) ' By=0.2778157677856496
    %tune B1 z0=100 z1=3000 initial=-0.6500 step=0.01 expr=Px1/Pz1 \
    %tolerance=0.000001
    %genericbend B fieldWidth=500 fieldHeight=500 fieldLength=500 \
    %ironWidth=800 ironHeight=800 ironLength=500 \
    %fieldMaterial=Vacuum place B z=2000 rename=B1 rotation=Y15 By=B1
    
    if sum(strcmp(char(B{1}(i)),defined))==0 % if not defined, define.
        Fw=240*1.0;
        Fh=70;
        Iw=10+Fw;
        Ih=30+Fh;
        
        split=1; %make 'split' parts of a single magnet
        defB=[ defB ...
            '\n param B' num2str(i) B{1}(i) '=' num2str(-Brho*B{3}(i)/B{4}(i),15) '\n'...
            ];
        if B{2}(i)>140 && B{2}(i)<1140
            
            defB=[ defB ...
                %  ' tune B' num2str(i) B{1}(i) ' z0=' num2str(B{2}(i)*1000-B{4}(i)*500) ' z1=' num2str(B{2}(i)*1000+B{4}(i)*500) ' '...
                %  ' initial=' num2str(Brho*B{3}(i)/B{4}(i),15) ' step=1e-13'...
                %  ' expr=(Px1/Pz1) tolerance=0.0000001 \n'...
                ];
        end
        defB=[ defB ...
            'genericbend ' B{1}(i) ' '...
            ' fieldLength=' num2str(B{4}(i)*1000/split,15) ... % mm
            ' fieldWidth=' num2str(Fw,15) ' '...
            ' fieldHeight=' num2str(Fh,15) ' '...
            ' ironLength=' num2str(B{4}(i)*1000/split,15) ... % mm
            ' ironWidth=' num2str(Iw,15) ' '...
            ' ironHeight=' num2str(Ih,15) ' '...
            ' fringe=0'...
            ' maxStep=' num2str(StepLength)  ...
            ... ' By=' num2str(Brho*B{3}(i)/B{4}(i),15)...' fieldColor=0.5,.1,.0 '...
            ' By=$B' num2str(i) B{1}(i) ...
            ' ironColor=1,.0,.0 kill=' num2str(killval) ' '...
            '\n'...
            ];
        % add to already defined
        defined{d}=char(B{1}(i));
        d=d+1;
        
    end
    
    plac{ord}=[
        'place ' B{1}(i) ' parent=  z=0+' num2str(0) '+0+' num2str(B{2}(i)*1000,15) ...
        ' x=' num2str(B{4}(i)/2*tan(B{3}(i)/2)*1000) ''... % move in by a sagitta amount
        ' coordinates=centerline'...' \n'...
        ' rotation=Y' num2str(B{3}(i)*180/pi/2,15)...' rename=' B{1}(i) ' \n'...
        ...' By=B' num2str(i) B{1}(i) ...
        '\n' ...
        ' cornerarc' ' z=0+' num2str(0) '+0+' num2str((B{2}(i)-B{4}(i)/2)*1000,15) ...
        ' angle=' num2str(B{3}(i)*180/pi,15)  ' centerRadius='  num2str((B{4}(i)/B{3}(i))*1000,15) '\n'...'#corner' ' parent=  z=' num2str((B{2}(i)+B{4}(i)/2)*1000,15) ... ' rotation=Y' num2str(B{3}(i)*180/pi/2,15) ' \n'...
        '\n']; %#ok<*SAGROW>
    ord=ord+1;
    
end

RF_on=0;
% plac RFC
%length[m]        voltage[MV]                lag          freq[MHz]             harmon
%0.501777                0.6           0.424389        476.0054439               1998
RF=['pillbox RFC maxGradient=' num2str(RF_on*0.6/0.501777,15) ...
    ' color=0.5,0.5,0 frequency=0.4760054439 innerLength=501.777 '...
    ' phaseAcc=0.424389 wallThick=100 innerRadius=1000\n'];
RF1pos=577.1465514*1000;
RF2pos=579.1935514*1000;
plac{ord}=['  # place RFC z=' num2str(RF1pos) '\n'];
plac{ord+1}=['  # place RFC z=' num2str(RF2pos) '\n'];

ippos=0;
Dumppos=beamstart-(+11+280);

positions=[0;Q{2}; S{2}; B{2}; (RF1pos-TotLength)/1000; (RF2pos-TotLength)/1000 ;ippos;(beamstart-TotLength)/1000;(TotLength/2)/1000;TotLength*2];
posString=[];
for i=1:length(positions)-1
    posString=[posString num2str(positions(i)*1000) ','];
end
posString=[posString num2str(positions(i)*1000)];

twissH=['profile ' ...
    'zloop=' num2str(beamstart*0) ':' num2str(TotLength*3) ':1000 '...
    ...' z=' posString ...
    ' file=twissLER particle=' particle ...
    ' coordinates=centerline\n']; %   ' coordinates=reference\n'];

%twissH=['profile zloop=' num2str(beamstart*0) ':' num2str(beamstart+TotLength*2.5) ':1000 file=twissHER particle=' particle ...
%      ' coordinates=centerline\n']; %   ' coordinates=reference\n'];


dump=['tubs DUMP innerRadius=0.01 outerRadius=300 length=10 kill=1 material=Cu color=0.2,0.2,0\n '...
    ];

% some hidrogen for interactions
ip=['material ebunch a=0.005 z=1 density=0.01\n'...
    'tubs ip innerRadius=0 outerRadius=300 length=1 kill=' num2str(killval) ' material=Vacuum color=1,.0,.1 \n '...
    ];
plac{ord+2}=['# place ip parent=  z=' num2str(ippos,15) ' \n' ];
plac{ord+3}=[...
    posbeam 'place BeamVis z=' num2str(beamstart,15) ' y=101 \n'
    ];
plac{ord+4}=['# place DUMP parent=  z=' num2str(Dumppos,15) ' \n'];
plac{ord+5}=[' cornerarc' ' z=0+' num2str(TotLength) '+0+' num2str(0,15) ...
    ' angle=' num2str(Xang/2*180/pi,15)  ' centerRadius='  num2str(1,15) '\n'...
    ];

%%%%%%                             RF1         RF2                   beam
[s,ind]=sort([0;Q{2}; S{2}; B{2}; RF1pos/1000; RF2pos/1000 ;637.2728766;beamstart/1000;Dumppos/1000;TotLength/1000]); % in m


HER=[head defQ defS defB RF ip dump];
HERPlace=[plac{ind}];

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       LER          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% LER (solenoids to implement)
folder='SBLERV12DATA';
P_0=4.179999969; % GeV/c
particle='e-';

Brho=3.3356409519815205*P_0; % electron

quadL=-1; % swich off LER quad.


killval=1;

beamstart=575.4563172*1000-140; % first QDI center

beamstart=0;

%beamstart=575.4563172*1000-140*3; % first QDI center
beamstart=(beamstart+TotLength);
%%%%%  %%%% %%%% both beams heading in the same direction!
elebeam=[ 'beam gaussian particle=' particle ' beamZ=' num2str(beamstart) ...
    ...' rotation=Y' num2str(1*Xang*180/pi) ... % if starting at IP.
    ...' rotation=Y' num2str(180) ... % if starting in straight.
    ' beamX=' num2str(-1*5.634489178e-05,15) ...
    ' beamXp=' num2str(+1*5.584726107e-08,15) ...
    '  nEvents=' num2str(Nevents) ' '... % nEvents=10
    ' sigmaX=' num2str(sX,15) ...
    ' sigmaY=' num2str(sY,15)...
    ' sigmaXp=' num2str(sXp,15)...
    ' sigmaYp=' num2str(sYp,15) ' '...
    '	meanMomentum=' num2str(P_0*1000,15) ' sigmaP=0'...
    ' meanT=0 sigmaT=0 \n'...
    'trackcolor reference=1,1,1 \n'...
    ];

headL=[...
    'reference particle=' particle ' referenceMomentum=' num2str(P_0*1000,15) ' noEloss=1'...
    ' beamXp=0 beamYp=0 '...
    ' beamZ=' num2str(TotLength) '+0 '...
    ' rotation=Y' num2str(0) ' '... % elements in the opposite direction
    ' tuneMomentum=' num2str(P_0,15) '\n'...
    ];


%headL=[...
%];
% x damping time= 0.26700741E-01 => 15% =~ 40ms = 40*10^6 ns

% get pos K and L for Quad, sext and bend
Qfile=fopen([folder '/quad']);
Sfile=fopen([folder '/sext']);
Bfile=fopen([folder '/bend']);
SOLfile=fopen([folder '/sol']);
% get pos K and L for Quad, sext and bend -- sbend and not rbend.
%Qfile=fopen('DIAMDATA/quad');
%Sfile=fopen('DIAMDATA/sext');
%Bfile=fopen('DIAMDATA/bend');

% name s K L
Q=textscan(Qfile,'%s %f %f %f','Headerlines',47);
B=textscan(Bfile,'%s %f %f %f','Headerlines',47);
S=textscan(Sfile,'%s %f %f %f','Headerlines',47);
SOL=textscan(SOLfile,'%s %f %f %f','Headerlines',47);

defQ= ['virtualdetector Det radius=10.0 length=' num2str(2*StepLength)...
    ' color=0,.2,.0 referenceParticle=1 format=rootExtended\n'];
defS='';
defSOL='';
defB='';
ord=1;
placL={};
defined={};
d=1;
% tilt reference by crossing angle.
placL{ord}=[' cornerarc' ' z=0+' num2str(TotLength) '+0+' num2str(0,15) ...
    ' angle=' num2str(Xang/2*180/pi,15)  ' centerRadius='  num2str(1,15) '\n'...
    ];
ord=ord+1;

for i=1:length(Q{1})
    name=(Q{1}(i));
    name{1}(1)=' ';
    name{1}(end)='p';
    Q{1}(i)=name;
end
for i=1:length(S{1})
    name=S{1}(i);
    name{1}(1)=' ';
    name{1}(end)='p';
    S{1}(i)=name;
end
for i=1:length(B{1})
    name=B{1}(i);
    name{1}(1)='';
    name{1}(end)='p';
    B{1}(i)=name;
end
for i=1:length(SOL{1})
    name=SOL{1}(i);
    name{1}(1)='';
    name{1}(end)='p';
    SOL{1}(i)=name;
end

% small FF quad
ffQuad{1}=char(Q{1}(1));%{'QD0AL'};
ffQuad{2}=char(Q{1}(2));%{'QD0BL'};
ffQuad{3}=char(Q{1}(3));%{'QD0AR'};
ffQuad{4}=char(Q{1}(4));
ffQuad{5}=char(Q{1}(end));%{'QD0BR'};
ffQuad{6}=char(Q{1}(end-1));%{'QDPMR'};
ffQuad{7}=char(Q{1}(end-2));%{'QDPML'};
ffQuad{8}=char(Q{1}(end-3));


for i=1:length(Q{1})
    %genericquad QDpi fieldLength=280 apertureRadius=101.5 ironRadius=381 \
    %	ironLength=250 ironColor=0.6,.0,.6 kill=' num2str(killval) ' gradient=$KQDpi
    %
    if sum(strcmp(char(Q{1}(i)),defined))==0;
        inrad=60;
        outrad=70;
        
        if sum(strcmp(char(Q{1}(i)),ffQuad))~=0
            inrad=12;
            outrad=15;
            tempSL=StepLength;
            StepLength=0.01;
            if strcmp(char(Q{1}(i)),ffQuad{1}) || strcmp(char(Q{1}(i)),ffQuad{5}) %% P.M.
                inrad=8;
                outrad=10;
                
            end
        end
        
        defQ=[ defQ ...
            ' genericquad ' Q{1}(i) ' '...
            ' fieldLength=' num2str(Q{4}(i)*1000,15) ... % mm
            ' apertureRadius=' num2str(inrad,15) ' '...
            ' ironRadius=' num2str(outrad,15) ' '...
            ' ironLength=' num2str(Q{4}(i)*1000,15) ' '...
            ' ironColor=0/255,255/255,191/255 kill=' num2str(killval) ' '...
            ' fringe=0'...
            ' gradient=' num2str(quadL*Q{3}(i)*Brho/Q{4}(i),15) ' '... K
            ' maxStep=' num2str(StepLength/2)  ...
            '\n']; %#ok<*AGROW>
        
        % add to already defined
        defined{d}=char(Q{1}(i));
        d=d+1;
    end
    
    StepLength=tempSL;
     
    %defQ=[defQ 'virtualdetector Det radius=400.0 color=0,.1,.1'];
    
    if (ord/2-floor(ord/2))==(det-1) % set 0 to place a detector in every quadrupole.
        placL{ord}=['place ' Q{1}(i) ' parent=  z=0+' num2str(TotLength) '+0+' num2str(Q{2}(i)*1000,15)...
            ' y=0'...
            ' coordinates=centerline'...
            '\n'...
            'place ' 'Det' ' parent=  z=0+' num2str(TotLength) '+0+' num2str((Q{2}(i))*1000)...
            ' rename=DETp# '...
            ' y=0'...
            ' coordinates=centerline '...
            '\n']; %#ok<*SAGROW>
    else
        placL{ord}=['place ' Q{1}(i) ' parent=  z=0+' num2str(TotLength) '+0+' num2str(Q{2}(i)*1000,15)...
            ' y=0'...
            ' coordinates=centerline'...
            '\n']; %#ok<*SAGROW>
    end
    %      if sum(strcmp(Q{1}(i),{'"QD0AL"';'"QD0BL"';'"QD0AR"';'"QD0BR"'}))>0
    %
    %      plac{ord}=['place ' Q{1}(i) ' parent=  z=' num2str(Q{2}(i)*1000,15)...
    %      ' coordinates=centerline'...
    %      ' rotation=Y' num2str(-0.066*180/pi,15) ' '...
    %      '\n'];
    %     end
    ord=ord+1;
end

for i=1:length(S{1})
    %multipole SF1 fieldLength=150 apertureRadius=101.5 ironRadius=381 \
    %	ironLength=140 ironColor=0,.6,.6 kill=' num2str(killval) ' sextupole=$KSF1
    
    if sum(strcmp(char(S{1}(i)),defined))==0;
        inrad=80;
        outrad=120;
        
        defS=[ defS ...
            ' multipole ' S{1}(i) ' '...
            ' fieldLength=' num2str(S{4}(i)*1000,15) ... % mm
            ' apertureRadius=' num2str(inrad,15) ' '...
            ' ironRadius=' num2str(outrad,15) ' '...
            ' ironLength=' num2str(S{4}(i)*1000,15) ' '...
            ' ironColor=173/255,255/255,47/255 kill=' num2str(killval) ' '...
            ' fringe=0'...
            ' maxStep=' num2str(StepLength/2)  ...
            ' sextupole=' num2str(sxt*S{3}(i)*Brho/S{4}(i),15) ' '... K2
            '\n'];
        % add to already defined
        defined{d}=char(S{1}(i));
        d=d+1;
    end
    
    
    placL{ord}=['place ' S{1}(i) ' parent=  z=0+' num2str(TotLength) '+0+' num2str(S{2}(i)*1000,15)...
        ' y=0'...
        ' coordinates=centerline '...
        '\n']; %#ok<*SAGROW>
    ord=ord+1;
    
end

for i=1:length(SOL{1})
    %multipole SF1 fieldLength=150 apertureRadius=101.5 ironRadius=381 \
    %	ironLength=140 ironColor=0,.6,.6 kill=' num2str(killval) ' sextupole=$KSF1
    
    if sum(strcmp(char(SOL{1}(i)),defined))==0;
        inrad=80;
        outrad=120;
        
        defSOL=[ defSOL ...
            ' \ncoil C' SOL{1}(i) ...
            ' innerRadius=' num2str(inrad) ...
            ' outerRadius=' num2str(outrad) ...
            ' length=' num2str(SOL{4}(i)) ...
            ' material=Cu' ...
            ' maxR=' num2str(outrad) '\n\n'...
            ' solenoid  ' SOL{1}(i) ' coil=C' SOL{1}(i) ' '...
            ' current=' num2str(0) ... %%% fake!!!!! not ok!!!
            ' color=0,1,1 kill=' num2str(killval) ...
            '\n\n'];
        % add to already defined
        defined{d}=char(SOL{1}(i));
        d=d+1;
    end
    
    
    placL{ord}=['place ' SOL{1}(i) ' parent=  z=0+' num2str(TotLength) '+0+' num2str(SOL{2}(i)*1000,15)...
        ' y=0'...
        ' coordinates=centerline '...
        '\n']; %#ok<*SAGROW>
    ord=ord+1;
    
end

for i=1:length(B{1})
    
    %genericbend BHer fieldWidth=1660 fieldHeight=200 fieldLength=5400 \
    %	ironColor=1,0,0 ironWidth=1828 ironHeight=1320 ironLength=5000 \
    % kill=' num2str(killval) ' By=0.2778157677856496
    if sum(strcmp(char(B{1}(i)),defined))==0 % if not defined, define.
        
        Fw=240;
        Fh=70;
        Iw=250;
        Ih=100;
        
        defB=[ defB ...
            '\n param Bp' num2str(i) B{1}(i) '=' num2str(Brho*B{3}(i)/B{4}(i),15) '\n'...
            ];
        if B{2}(i)<400 || B{2}(i)>800
            defB=[ defB ...
                %  ' tune Bp' num2str(i) B{1}(i) ' z0=' num2str(B{2}(i)*1000-B{4}(i)*500) ' z1=' num2str(B{2}(i)*1000+B{4}(i)*500) ' '...
                %  ' initial=' num2str(Brho*B{3}(i)/B{4}(i),15) ' step=1e-13'...
                %  ' expr=(Px1/Pz1) tolerance=0.0000001 \n'...
                ];
        end
        defB=[ defB ...
            ' genericbend ' B{1}(i) ' '...
            ' fieldLength=' num2str(B{4}(i)*1000,15) ... % mm
            ' fieldWidth=' num2str(Fw,15) ' '...
            ' fieldHeight=' num2str(Fh,15) ' '...
            ' ironLength=' num2str(B{4}(i)*1000,15) ... % mm
            ' ironWidth=' num2str(Iw,15) ' '...
            ' ironHeight=' num2str(Ih,15) ' '...
            ' ironLength=' num2str(B{4}(i)*1000,15) ' '...
            ' fringe=0'...
            ' maxStep=' num2str(StepLength)  ...
            ...' By=' num2str(Brho*B{3}(i)/B{4}(i),15)...' fieldColor=0.5,.1,.0 '...
            ' By=$Bp' num2str(i) B{1}(i) ...
            ' ironColor=1,.0,.0 kill=' num2str(killval) ' '...
            '\n'...
            ];
        % add to already defined
        defined{d}=char(B{1}(i));
        d=d+1;
        
    end
    
    
    %'# cornerarc' ' parent=  z=' num2str((B{2}(i))*1000,15) ...    ' angle=' num2str(B{3}(i)*180/pi) ' centerRadius='  num2str((B{4}(i)/B{3}(i))*1000,15) '\n'...'# cornerarc' ' parent=  z=' num2str((B{2}(i)-B{4}(i)/2)*1000,15) ...    ' angle=' num2str(B{3}(i)*180/pi/2,15) ' centerRadius='  num2str((B{4}(i)/B{3}(i))*1000,15) '\n'...'#corner' ' parent=  z=' num2str((B{2}(i)-B{4}(i)/2)*1000,15) ...' rotation=Y' num2str(B{3}(i)*180/pi/2,15) ' \n'...
    placL{ord}=[
        'place ' B{1}(i) ' parent=  z=0+' num2str(TotLength) '+0+' num2str(B{2}(i)*1000,15) ...
        ' x=' num2str(B{4}(i)/2*tan(B{3}(i)/2)*1000) ''...
        ' coordinates=centerline'...' \n'...
        ' rotation=Y' num2str(B{3}(i)*180/pi/2,15) ...' rename=' B{1}(i) ' \n'...
        '\n' ...
        ' cornerarc' ' z=0+' num2str(TotLength) '+0+' num2str((B{2}(i)-B{4}(i)/2)*1000,15) ...
        ' angle=' num2str(B{3}(i)*180/pi,15)  ' centerRadius='  num2str((B{4}(i)/B{3}(i))*1000,15) '\n'...'#corner' ' parent=  z=' num2str((B{2}(i)+B{4}(i)/2)*1000,15) ... ' rotation=Y' num2str(B{3}(i)*180/pi/2,15) ' \n'...
        '\n']; %#ok<*SAGROW>
    ord=ord+1;
    
end

RF_on=0;

% plac RFC
%length[m]        voltage[MV]                lag          freq[MHz]             harmon
%0.501777                0.6           0.424389        476.0054439               1998
RFL=['pillbox RFC maxGradient=' num2str(RF_on*0.6/0.501777,15) ...
    ' color=0.5,0.5,0 frequency=0.4760054439 innerLength=501.777 '...
    ' phaseAcc=0.424389 wallThick=100 innerRadius=1000\n'];
RF1pos=TotLength+577.1465514*1000;
RF2pos=TotLength+579.1935514*1000;

placL{ord}=['  # place RFC z=' num2str(RF1pos) '\n'];
placL{ord+1}=['  # place RFC z=' num2str(RF2pos) '\n'];

%rfc, at = 11.2402595;
%rfc, at = 13.2872595;
positions=[0;Q{2}; S{2}; SOL{2}; B{2}; (RF1pos-TotLength)/1000; (RF2pos-TotLength)/1000 ;ippos;(beamstart-TotLength)/1000;(TotLength/2)/1000;TotLength*2];
posString=[];
for i=1:length(positions)-1
    posString=[posString num2str(positions(i)*1000+TotLength) ','];
end
posString=[posString num2str(positions(i)*1000+TotLength)];

twissL=['profile ' ...
    'zloop=' num2str(beamstart*0) ':' num2str(TotLength*3) ':1000 '...
    ...' z=' posString ...
    ' file=twissLER particle=' particle ...
    ' coordinates=centerline\n']; %   ' coordinates=reference\n'];

% keep only part of sequence.
% ind=ind(end/4:end/4*2);

dumpL=['tubs DUMP innerRadius=0 outerRadius=300 length=10 kill=1 material=Cu color=0.2,0.2,0\n '...
    ];

% some hidrogen for interactions
ipL=['material ebunch a=0.005 z=1 density=0.01\n'...
    'tubs ip innerRadius=0 outerRadius=300 length=1 kill=' num2str(killval) ' material=Vacuum color=1,.0,.1 \n '...
    ];
ippos=0;
Dumppos=beamstart+(-11-280);
placL{ord+2}=['# place ip parent=  z=' num2str(ippos,15) ' \n' ];
placL{ord+3}=[elebeam 'place BeamVis z=' num2str(beamstart,15) ' y=-101 \n'];
placL{ord+4}=['# place DUMP parent=  z=' num2str(Dumppos,15) ' \n'];

placL{ord+5}=[' cornerarc' ' z=0+' num2str(TotLength*2) '+0+' num2str(0,15) ...
    ' angle=' num2str(-Xang/2*180/pi,15)  ' centerRadius='  num2str(1,15) '\n'...
    ];


% first 0 is for cornerarc crossing angle rotation.
[s,ind]=sort([0;Q{2}; S{2}; SOL{2}; B{2}; (RF1pos-TotLength)/1000; (RF2pos-TotLength)/1000 ;ippos;(beamstart-TotLength)/1000;(TotLength/2)/1000;TotLength*2]); % in m



outtext=char(cell2mat([...
    bigmondo ...
    ... definitions
    HER ...  HER
    defQ defS defB defSOL ...  LER
    RF ip...
    HERPlace ...  Place HER magnets
    twissH ...
    headL  dumpL...
    placL{ind} twissL...  Place LER magnets
    ' particlecolor e-=0,0,1 e+=1,0,0 \n'...  electron is blue and positron is red.
    ...'g4ui when=4 ''/vis/open OGL''\n'...
    'g4ui when=4 ''/vis/viewer/set/viewpointThetaPhi 45 45''\n'...
    'g4ui when=4 ''/vis/viewer/set/style wireframe''\n'...
    'g4ui when=4 ''/run/beamOn ' num2str(Nevents*2) ' ''\n'...
    ... 'g4ui when=4 ''vis/ogl/printEPS''\n'...
    ]));
% lattice 'g4ui when=4 ''/vis/viewer/set/style wireframe''\n'...    'g4ui ''/run/beamOn 20''\n'...
out=fopen('G4blseqconv','w+');
fprintf(out,outtext);

fclose('all');
clear defQ
clear defS
clear defB
clear plac ind ord
clear defined
disp('--- DONE ---')