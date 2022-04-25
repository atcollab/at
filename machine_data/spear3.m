function varargout = spear3(varargin)
%SPEAR3 Load the SPEAR3 lattice structure

L0  = 2.341440122400003e+002;	% design length [m] 
C0  = 299792458; 				% speed of light [m/s]
H   = 372;                      % Cavity harmonic number

CAV	= struct('FamName', 'RF' , 'Energy', 3e9, ...
        'Length', 0 , 'Voltage', 3.2e+6 , 'Frequency', H*C0/L0, ...
        'HarmNumber', H , 'PassMethod', 'RFCavityPass');

    
COR = struct('FamName', 'COR' ,...
        'Length', 0.15, 'KickAngle', [0 0], 'PassMethod', 'CorrectorPass'); 
    
BPM.FamName = 'BPM';
BPM.Length = 0;
BPM.PassMethod = 'IdentityPass';

AP.FamName = 'AP';
AP.Length = 0;
AP.Limits = [-0.1, 0.1, -0.1, 0.1];
AP.PassMethod = 'AperturePass';

% ===================== Injection Kickers and Drifts ============================================ 

INJ.FamName = 'SEPTUM';
INJ.MemberOf = 'Injection';
INJ.Length = 0;
INJ.PassMethod = 'IdentityPass'; 

K1 = struct('FamName', 'KICKER' , 'MemberOf', {{'Injection'}}, 'Tag', 'K1', ...
        'Length', 1.2 , 'KickAngle', [0 0], 'PassMethod', 'CorrectorPass'); 
    

K2 = K1; 
K2.Length = 0.6;
K2.Tag = 'K2';


K3 = K1;
K3.Tag = 'K3';


DI1 = atelem('drift', 'FamName', 'DI1'  ,'Length', 0.9235741);
DI2 = atelem('drift', 'FamName', 'DI2'  ,'Length', 0.6882939);
DI3 = atelem('drift', 'FamName', 'DI3'  ,'Length', 0.6834939);
DI4 = atelem('drift', 'FamName', 'DI4'  ,'Length', 0.1224401);
DI5 = atelem('drift', 'FamName', 'DI5'  ,'Length', 1.24013);
DI6 = atelem('drift', 'FamName', 'DI6'  ,'Length', 0.165804);



% ================  Standard Cell Drifts ===================================

DC2 = atelem('drift', 'FamName', 'DC2', 'Length',  0.097500);
DC5 = atelem('drift', 'FamName', 'DC5', 'Length',  0.200986);

% ================= Standard Cell BPM Drifts ===============================

DC1A    =    atelem('drift', 'FamName', 'DC1A', 'Length', 1.405934);
DC1B    =    atelem('drift', 'FamName', 'DC1B', 'Length', 0.12404125);
DC3A    =    atelem('drift', 'FamName', 'DC3A', 'Length', 0.05322065);
DC3B    =    atelem('drift', 'FamName', 'DC3B', 'Length', 0.16368247);

DC4A    =    atelem('drift', 'FamName', 'DC4A', 'Length', 0.15921467);
DC4B    =    atelem('drift', 'FamName', 'DC4B', 'Length', 0.044418);
DC6A    =    atelem('drift', 'FamName', 'DC6A', 'Length', 0.110646);
DC6B    =    atelem('drift', 'FamName', 'DC6B', 'Length', 0.06316585);  %0.069354 corrected to make path length consistent

% ================= Standard Cell Corrector Magnet Drifts ==================
DC2A    =    atelem('drift', 'FamName', 'DC2A', 'Length', 0.11576525);
DC2B    =    atelem('drift', 'FamName', 'DC2B', 'Length', 0.11581045);
DC2C    =    atelem('drift', 'FamName', 'DC2C', 'Length', 0.10210045);
DC2D    =    atelem('drift', 'FamName', 'DC2D', 'Length', 0.12947525);

DC5A    =    atelem('drift', 'FamName', 'DC5A', 'Length', 0.09058);
DC5B    =    atelem('drift', 'FamName', 'DC5B', 'Length', 0.36139);
DC5C    =    atelem('drift', 'FamName', 'DC5C', 'Length', 0.09584);
DC5D    =    atelem('drift', 'FamName', 'DC5D', 'Length', 0.35613);


% ================ Bending Magnets ======================================

BEND	=	atelem('bend', 'FamName', 'BND','Length', 1.5048,  ...
            'BendingAngle', 0.18479957, 'EntranceAngle', 0.18479957/2,...
            'ExitAngle', 0.18479957/2, 'K', -0.31537858);
 
BDM	    =	atelem('bend'  , 'FamName', 'B34', 'Length', 1.14329,  ...
            'BendingAngle', 0.138599675894, 'EntranceAngle', 0.138599675894/2, ...
            'ExitAngle', 0.138599675894/2, 'K', -0.31537858);
        
% ================ Standard Cell Quadrupoles  ===========================

QF		=	atelem('quadrupole', 'FamName', 'QF', 'Length', 0.3533895,  'K', 1.768672904054 );
QD		=	atelem('quadrupole', 'FamName', 'QD', 'Length', 0.1634591,  'K', -1.542474230359 );
QFC	    =	atelem('quadrupole', 'FamName', 'QFC', 'Length', 0.5123803, 'K', 1.748640831069);

% ================ Matching Cell Quadrupoles  ===========================

QDX 	=	atelem('quadrupole', 'FamName', 'QDX', 'Length', 0.3533895,  'K', -1.386467245226 );
QDY     =   atelem('quadrupole', 'FamName', 'QDY', 'Length', 0.3533895,  'K', -0.460640930646 );
QDZ     =   atelem('quadrupole', 'FamName', 'QDZ', 'Length', 0.3533895,  'K', -0.878223937747 );
QFX 	=	atelem('quadrupole', 'FamName', 'QFX', 'Length', 0.6105311,  'K', 1.573196278394 );
QFY 	=	atelem('quadrupole', 'FamName', 'QFY', 'Length', 0.5123803,  'K', 1.481493709831 );
QFZ 	=	atelem('quadrupole', 'FamName', 'QFZ', 'Length', 0.3533895,  'K', 1.427902006984 );


% ================ Sextupoles ============================================

SF = atelem('sextupole', 'FamName', 'SF' , 'Length' , 0.21);
SF.PolynomB(3) = 32.0477093/2;
SF.MaxOrder = 2;

SFM	 =	SF; % SFM is the same length as SF
SFM.FamName = 'SFM';
SFM.PolynomB(3) = 7.5;

SDM	=	SF; % SFM is the same length as SF, SFM
SDM.FamName = 'SDM';
SDM.PolynomB(3) = -8.5;

SD = atelem('sextupole', 'FamName', 'SD' , 'Length' , 0.25);
SD.PolynomB(3) = -38.80153/2;
SD.MaxOrder = 2;

% ============== Matching Cell Drifts without correctors or BPMs
% NOTE: BPMS and correctors are not symmetric in MCA, MCB

DM1	    =	atelem('drift', 'FamName','DM1' ,'Length', 3.81);
DM2	    =	atelem('drift', 'FamName','DM2' ,'Length', 0.0975);
DM3	    =	atelem('drift', 'FamName','DM3' ,'Length', 0.275);
DM4	    =	atelem('drift', 'FamName','DM4' ,'Length', 0.21584572);
DM5	    =	atelem('drift', 'FamName','DM5' ,'Length', 0.250);
DM6	    =	atelem('drift', 'FamName','DM6' ,'Length', 0.49068463);
DM7	    =	atelem('drift', 'FamName','DM7' ,'Length', 0.17380985);
DM8 	=	atelem('drift', 'FamName','DM8' ,'Length', 0.500);
DM9	    =	atelem('drift', 'FamName','DM9' ,'Length', 0.105);
DM10    =	atelem('drift', 'FamName','DM10','Length', 3.2765714);

%Matching Cell A BPM Drifts
DA1A	=	atelem('drift', 'FamName','DA1A' , 'Length' ,3.6792386);
DA1B	=	atelem('drift', 'FamName','DA1B' , 'Length' ,0.12406665);

DA3A	=	atelem('drift', 'FamName','DA3A' , 'Length' ,0.20889925);
DA3B	=	atelem('drift', 'FamName','DA3B' , 'Length' ,0.05414045);

DA5A	=	atelem('drift', 'FamName','DA5A' , 'Length' ,0.11397747);
DA5B	=	atelem('drift', 'FamName','DA5B' , 'Length' ,0.108563);
DA5C	=	atelem('drift', 'FamName','DA5C' , 'Length' ,0.051845);   
DA5D	=	atelem('drift', 'FamName','DA5D' , 'Length' ,0.17069547);   

DA7A	=	atelem('drift', 'FamName','DA7A' , 'Length' ,0.1106966);
DA7B	=	atelem('drift', 'FamName','DA7B' , 'Length' ,0.06311325);

DA8A	=	atelem('drift', 'FamName','DA8A' , 'Length' ,0.33735947);
DA8B	=	atelem('drift', 'FamName','DA8B' , 'Length' ,0.12848625);

DA10A   =	atelem('drift', 'FamName','DA10A', 'Length' ,0.12393965);
DA10B   =	atelem('drift', 'FamName','DA10B', 'Length' ,3.145937 );

%Matching Cell A Corrector Drifts
DA2A	=	atelem('drift', 'FamName','DA2A' , 'Length' ,0.11530525);
DA2B	=	atelem('drift', 'FamName','DA2B' , 'Length' ,0.11773445);

DA6A	=	atelem('drift', 'FamName','DA6A' , 'Length' ,0.1266);
%DA6B	=	atelem('drift', 'FamName','DA6B' , 'Length' ,0.90476852);   %0.90477 corrected to make path length consistent
DA6B	=	atelem('drift', 'FamName','DA6B' , 'Length' ,0.90476828);   %0.90477 corrected to make path length consistent with MAD (234.14401272)
DA6C	=	atelem('drift', 'FamName','DA6C' , 'Length' ,0.0960);
DA6D	=	atelem('drift', 'FamName','DA6D' , 'Length' ,0.93537);

DA9A	=	atelem('drift', 'FamName','DA9A' , 'Length' ,0.10930525);
DA9B	=	atelem('drift', 'FamName','DA9B' , 'Length' ,0.13730525);

%Matching Cell B  BPM Drifts 
DB1A	=	atelem('drift', 'FamName','DB1A' , 'Length' ,3.747082 );
DB1B	=	atelem('drift', 'FamName','DB1B' , 'Length' ,0.05622325 );

DB3A	=	atelem('drift', 'FamName','DB3A' , 'Length' ,0.13222685);
DB3B	=	atelem('drift', 'FamName','DB3B' , 'Length' ,0.13081285);

DB5A	=	atelem('drift', 'FamName','DB5A' , 'Length' ,0.17069547);
DB5B	=	atelem('drift', 'FamName','DB5B' , 'Length' ,0.051845 );
DB5C	=	atelem('drift', 'FamName','DB5C' , 'Length' ,0.1085632);
DB5D	=	atelem('drift', 'FamName','DB5D' , 'Length' ,0.11397727);

DB7A	=	atelem('drift', 'FamName','DB7A' , 'Length' ,0.06311305);
DB7B	=	atelem('drift', 'FamName','DB7B' , 'Length' ,0.1106968);

DB8A	=	atelem('drift', 'FamName','DB8A' , 'Length' ,0.32725027);
DB8B	=	atelem('drift', 'FamName','DB8B' , 'Length' ,0.13859545);

DB10A	=	atelem('drift', 'FamName','DB10A', 'Length' ,0.12404125 );
DB10B	=	atelem('drift', 'FamName','DB10B', 'Length' ,3.1458354);

%Matching Cell B Corrector Drifts
DB2A	=	atelem('drift', 'FamName','DB2A' , 'Length' ,0.115805250);
DB2B	=	atelem('drift', 'FamName','DB2B' , 'Length' ,0.117234450);

DB6A	=	atelem('drift', 'FamName','DB6A' , 'Length' ,0.93737);
DB6B	=	atelem('drift', 'FamName','DB6B' , 'Length' ,0.09399852);   %0.094 corrected to make path length consistent
DB6C	=	atelem('drift', 'FamName','DB6C' , 'Length' ,0.90437);
DB6D	=	atelem('drift', 'FamName','DB6D' , 'Length' ,0.1270);

DB9A	=	atelem('drift', 'FamName','DB9A' , 'Length' ,0.12330525);
DB9B	=	atelem('drift', 'FamName','DB9B' , 'Length' ,0.12330525);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%Standard Cell: (Note QFC not split)
HCEL1 =	{DC1A BPM DC1B QF DC2A COR DC2B QD DC3A BPM DC3B BEND DC4A BPM DC4B SD...
DC5A COR DC5B SF DC6A BPM DC6B QFC};

HCEL2 =	{DC6B DC6A SF DC5C COR DC5D...
SD DC4B BPM DC4A BEND DC3B DC3A QD DC2C COR DC2D QF DC1B BPM DC1A};

ACEL	=	[{AP} HCEL1 HCEL2];

%Cell 2: (K1 magnet, Note QFC not split)
K1CEL2 =	{DC6B DC6A SF DC5C COR DC5D...
SD DC4B BPM DC4A BEND DC3B DC3A QD DC2C COR DC2D QF DC1B BPM DI1};
CEL2	=	[{AP} HCEL1 K1CEL2];

%Cell 3: (K1 & K2 magnets, Note QFC not split)
K1CEL3 =	{K1 DI2 BPM DC1B QF DC2A COR DC2B QD DC3A BPM DC3B BEND DC4A BPM DC4B SD...
DC5A COR DC5B SF DC6A BPM DC6B QFC};

K2CEL3 =	{DC6B DC6A SF DC5C COR DC5D...
SD DC4B BPM DC4A BEND DC3B DC3A QD DC2C COR DC2D QF DC1B BPM DI3 K2 DI4};

CEL3	=	[{AP} K1CEL3 K2CEL3];

%Cell 4: (Septum & K3 magnets, Note QFC not split)
SEPCEL4 =	{DI5 INJ DI6 BPM DC1B QF DC2A COR DC2B QD DC3A BPM DC3B BEND DC4A BPM DC4B SD...
DC5A COR DC5B SF DC6A BPM DC6B QFC};

K3CEL4 =	{DC6B DC6A SF DC5C COR DC5D...
SD DC4B BPM DC4A BEND DC3B DC3A QD DC2C COR DC2D QF DC1B BPM DI2 K3};

CEL4	=	[{AP} SEPCEL4 K3CEL4];

%Cell 5: (K5 magnets, Note QFC not split)
K3CEL5 =	{DI1 BPM DC1B QF DC2A COR DC2B QD DC3A BPM DC3B BEND DC4A BPM DC4B SD...
DC5A COR DC5B SF DC6A BPM DC6B QFC};
CEL5	=	[{AP} K3CEL5 HCEL2];

%Matching Cell A (South East, North West)
MCA={DA1A BPM DA1B QDX DA2A COR DA2B ...
QFX DA3A BPM DA3B QDY DM4 BDM DA5A BPM DA5B SDM DA6A COR DA6B SFM DA7A BPM DA7B QFY ...
DM7 SFM DA6C COR DA6D SDM DA5C BPM DA5D BDM...
DA8A BPM DA8B QDZ DA9A COR DA9B QFZ DA10A BPM DA10B};

%Matching Cell B (North East, South West)
MCB={DB1A BPM DB1B QDX DB2A COR DB2B ...
QFX DB3A BPM DB3B QDY DM4 BDM DB5A BPM DB5B SDM DB6A COR DB6B SFM DM7 QFY DB7A BPM...
DB7B SFM DB6C COR DB6D SDM DB5C BPM DB5D BDM...
DB8A BPM DB8B QDZ DB9A COR DB9B QFZ DB10A BPM DB10B};

% Begin Lattice
NORTH   =	[CEL2 CEL3 CEL4 CEL5 ACEL ACEL ACEL];
SOUTH	=	[ACEL ACEL ACEL ACEL ACEL ACEL ACEL];
RING    =   [{CAV} MCA NORTH flip(MCB) MCA SOUTH flip(MCB)];

if nargout > 0
    varargout{1} = RING(:);
else % If no output arguments - greate global variable THERING
    global THERING
    if ~isempty(THERING)
        warning('Global variable THERING was overridden');
    end
    THERING = RING(:);
    evalin('base','global THERING');
end
