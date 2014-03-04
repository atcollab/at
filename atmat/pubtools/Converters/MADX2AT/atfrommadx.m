function atfrommadx(seqfilemadX,E0,outfilename)
%function atfrommadx(seqfilemadX,E0,outfilename)
% tansform madX sequence file (savesequence) file into AT lattice structure.
%
% This procedure reads a saved lattice (sequence in madx) in madX 
% and converts it to an AT lattice
%
% (madx comands to save the sequences : 
% 
%  _______ MADX code _________
%  use,period=sequencename1; 
%  use,period=sequencename2; 
%  use,period=sequencename2; 
%  SAVE,FILE='seqfilemadX.seq';
%  ___________________________
%
%  seqfilemadX.seq will contain sequencename1 sequencename2 sequencename3 
%  in the correct format in a single file
% 
% )
% 
% The routine outputs a Matlab macro with all the AT defitions and variables as
% in the madX file  
%
% The order of the declarations is the same in the two files.
% declarations that contain other variables are moved to the end. (this may not be enough)
%
%
% Works also with single madX files not containing comands, only
% definitions.
%
% parameters:
%    - seqfilemadX=name of the mad8 lattice file
%    - E0  = design energy
%    - outfilename (default: seqfilemadX_AT_LATTICE.mat)
%    
% default pass methods:
%          quadrupoles : StrMPoleSymplectic4Pass
%          dipole : BndMPoleSymplectic4Pass
%          multipole : StrMPoleSymplectic4Pass
%          sextupole : StrMPoleSymplectic4Pass
%          thinmultipole : ThinMPolePass
%          correctors : ThinMPolePass
%          cavity : DriftPass
%

%% changes history
% created 7-sep-2012 by S.M.Liuzzo @ ESRF
%
% updated 12-sep-2012 Cavity DriftPass (not IdentityPass)
% updated 13-sep-2012 Capital letter names
%                     CorrectorPass
%                     _red reduced versions
% updated 14-sep-2012 multipoles (1/n!) factor from madx to AT
% updated 21-dec-2012 rbend-sbend length conversion corrected
% updated 05-mar-2013 lattice output no fringes and with fringes (_FF).
%                     
% updated 20-mar-2013 removed output no fringes and with fringes (_FF).
%                     removed output reduced (_red).
%                     use atconstructors.
%                     call macro directly
%                     try-catch macro running
%                     tempname file and fulfile.
                  

%% initial warnings
disp('important notice about conversion:')
disp('');
disp(['1) THE MADX FILE MUST BE OUTPUT OF save, file=seqfilemadX;'...
    ' THIS GARANTES AN APPROPRIATE FILE FORMAT']);
disp('');
disp(['2) THE MADX PROGRAM allows to use variable not previously defined.'...
    ' This is not known to AT!.'...
    ' If the program fails but generates a ..._macro.m file, '...
    ' please edit this file reordering the declarations and run it.']);
disp('');
disp(['3) If periodname is not specified' ...
    ' It is assumed the name of the file (up to the .) ' ...
    'to be the sequence name that you want to use']);

%% open madX sequence file
sX=fopen(seqfilemadX,'r');

% the : detect definitions. every : is preceded by a name
% of an element.
% between two : all the parameter of an elemetn or a line ar found.

%% get madX file in a cell array with a new elment at every space comma or new line.
SXCELL=textscan(sX,'%s','delimiter',';');

SXSTRING=SXCELL{1}; % still a cell array but every space or new line now is stored as a new cell.

% scroll and reshape to divide atributes
tmp={};

for i=1:length(SXSTRING)
    SXSTRING{i}=strrep(SXSTRING{i},' ',''); % no spaces
    SXSTRING{i}=strrep(SXSTRING{i},':=','=');
    SXSTRING{i}=strrep(SXSTRING{i},';','');
    
    c=[1 sort([strfind(SXSTRING{i},',') strfind(SXSTRING{i},':')...
                     strfind(SXSTRING{i},'=')]) length(SXSTRING{i})];
    
    for jc=1:length(c)-1
        if jc==1
            tmp=[tmp SXSTRING{i}(c(jc):c(jc+1))];
        else
            tmp=[tmp SXSTRING{i}(c(jc)+1:c(jc+1))];
        end
    end
    
end
SXSTRING=tmp;


%% open .m file to output matlab translated code

filemacroname=[tempname '.m'];

mafileout=fopen(filemacroname,'w+');
mst=['%% this is a macro that converts to AT the madX lattice: '...
                     seqfilemadX '\n%%\n%% Created: ' datestr(now)...
                     '\n%%\n%%\n\nglobal GLOBVAL;\nGLOBVAL.E0=' num2str(E0)...
                     ';\n\n\n'];
def=['\n\n%%%% DEFINITIONS \n\n']; %#ok<*NBRAK>
lines=['\n\n%%%% LINES \n\n'];
var=['%%%% VARIABLES \n\n sxt_on=1;'];
formulas=['\n\n%%%% RELATIONS \n\n'];

%% convert to a matlab macro
j=1; % used in line block counter (element atributes counter)

elemcount=0;
i=1; % skip header in mad8 file   (element counter)
while i<length(SXSTRING)-2
    
    if SXSTRING{i}(end)==':' % new element or line
        def=[ def ...
            ]; %#ok<*AGROW> % end line and go to newline add name
        SXSTRING{i}(1:end-1)=upper(SXSTRING{i}(1:end-1));
        
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'MARKER','MARKER,');
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'MARKER,,','MARKER,');
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'KICKER','KICKER,');
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'KICKER,,','KICKER,');
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'MONITOR','MONITOR,');
        SXSTRING{i+j}=strrep(SXSTRING{i+j},'MONITOR,,','MONITOR,');
        
        ElementType=SXSTRING{i+j}(1:end-1);
        %  THERING=[THERING ' ' SXSTRING{i+1}];
        %j=2;
        nwel=SXSTRING{i+j}(end);
      
        switch ElementType
            case {'quadrupole','QUADRUPOLE', 'QUADRUPO'}
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atquadrupole('''...
                    SXSTRING{i}(1:end-1)...
                    ''',0,0,''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                   
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
            case {'sextupole','SEXTUPOLE'}
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atsextupole('''...
                    SXSTRING{i}(1:end-1)...
                    ''',0,0,''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')=' SXSTRING{i}(1:end-1) '.(''PolynomB'')*1/2' ';\n' ...
                ];
                
            case {'rbend','RBEND'}
                
                def=[def '\n'];
               
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atrbend('''...
                    SXSTRING{i}(1:end-1)...
                    ''',0,0,0,''BndMPoleSymplectic4Pass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                % bendings are  sector by default. change to rectangular
                def=[ def ...
                 SXSTRING{i}(1:end-1) '.(''EntranceAngle'')=' SXSTRING{i}(1:end-1) '.(''EntranceAngle'')+' SXSTRING{i}(1:end-1) '.(''BendingAngle'')/2; \n'...
                 SXSTRING{i}(1:end-1) '.(''ExitAngle'')=' SXSTRING{i}(1:end-1) '.(''ExitAngle'')+' SXSTRING{i}(1:end-1) '.(''BendingAngle'')/2; \n'...
                   SXSTRING{i}(1:end-1) '.(''Length'')=' SXSTRING{i}(1:end-1)...
                     '.(''Length'')*(' SXSTRING{i}(1:end-1)...
                     '.(''BendingAngle'')/2)/sin(' SXSTRING{i}(1:end-1)...
                     '.(''BendingAngle'')/2); \n'...
                  SXSTRING{i}(1:end-1) '.(''MaxOrder'')=length(' SXSTRING{i}(1:end-1) '.(''PolynomB''))-1; \n'];
              
            case {'sbend','SBEND'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atsbend('...
                    '''' SXSTRING{i}(1:end-1)...
                    ''',0,0,0,''BndMPoleSymplectic4Pass'');\n'...
                    ];
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
           def=[ def SXSTRING{i}(1:end-1) '.(''MaxOrder'')=length(' SXSTRING{i}(1:end-1) '.(''PolynomB''))-1; \n'];
              
            case {'DRIFT','drift'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atdrift('...
                    '''' SXSTRING{i}(1:end-1)...
                    ''',0,''DriftPass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'RFCAVITY','rfcavity'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atrfcavity(''' SXSTRING{i}(1:end-1)...
                    ''',0,0,0,0,' num2str(E0) ...
                    ',''DriftPass'');\n'...
                    ];
                
                elemcount=elemcount+1;
                
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'MULTIPOLE','multipole'}
                % multipoles should be StrMPoleSymplectic4Pass with short length
                % to be compatible with MADX
                
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '=atthinmultipole('...
                    '''' SXSTRING{i}(1:end-1) ''''...
                    ',[0 0 0 0],[0 0 0 0],'...
                    '''ThinMPolePass'');\n'...
                     SXSTRING{i}(1:end-1) '.(''Length'')=0; \n'...
                    ];
                
                
                elemcount=elemcount+1;
                % MADX-->   ocf0: multipole,knl:={ 0, 0, 0,kocf0 };
                while nwel~=':' % loops atributes of this element definition
                    multipoles=[];
                    if strcmp(SXSTRING{i+j+1}(1),'{') % if open paerntesis found 
                        multipoles=[multipoles '[' SXSTRING{i+j+1}(2:end)];
                        k=2;
                        while ~strcmp(SXSTRING{i+j+k}(end),'}') && k<10 % look for closed parentesis
                            multipoles=[multipoles SXSTRING{i+j+k}];
                            k=k+1;
                        end
                        multipoles=[multipoles SXSTRING{i+j+k}(1:(end-1)) ']'];
                        
                    end
                    
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},multipoles);
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''Class'')=''Multipole''; \n'...% max order size polynomb -1
                    SXSTRING{i}(1:end-1) '.(''MaxOrder'')=numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')' ')-1; \n'...% max order size polynomb -1
                    'expansionCoefFixA=1./factorial([1: numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomA''))]-1); \n'...
                    'expansionCoefFixB=1./factorial([1: numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB''))]-1); \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')=(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')' ').*expansionCoefFixB; \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomA'')=(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomA'')' ').*expansionCoefFixA; \n'...
                    ];
                
            case {'OCTUPOLE','octupole'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '=atmultipole('...
                    '''' SXSTRING{i}(1:end-1) ''''...
                    ',0,[0 0 0 0],[0 0 0 0],'...
                    '''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
                % fix madX-AT multipole coefficents
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''Class'')=''Octupole''; \n'...% max order size polynomb -1
                    SXSTRING{i}(1:end-1) '.(''MaxOrder'')=numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')' ')-1; \n'...% max order size polynomb -1
                    'expansionCoefFixA=1./factorial([1: numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomA''))]-1); \n'...
                    'expansionCoefFixB=1./factorial([1: numel(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomB''))]-1); \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')=(' ....
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')' ').*expansionCoefFixB; \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomA'')=(' ...
                    SXSTRING{i}(1:end-1) '.(''PolynomA'')' ').*expansionCoefFixA; \n'...
                    ];
                
            case {'SOLENOID','solenoid'}
                def=[def '\n'];
                
                 def=[ def ...
                    SXSTRING{i}(1:end-1) '=atsolenoid('...
                    '''' SXSTRING{i}(1:end-1) ''''...
                    ',0,0,'...
                    '''IdentityPass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'MARKER','marker','marke'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '=atmarker('...
                    '''' SXSTRING{i}(1:end-1) ''');\n'...
                     SXSTRING{i}(1:end-1) '.(''Length'')=0; \n'...
                     SXSTRING{i}(1:end-1) '.(''Class'')=''Marker''; \n'...
                    ];
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'MONITOR','monitor','monito'}
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''FamName'')=''' SXSTRING{i}(1:end-1) ''';\n' ...
                    SXSTRING{i}(1:end-1) '.(''Class'')=''Monitor''; \n'...
                    SXSTRING{i}(1:end-1) '.(''BetaCode'')=''PU''; \n'...
                    SXSTRING{i}(1:end-1) '.(''PassMethod'')=''IdentityPass''; \n'...
                    SXSTRING{i}(1:end-1) '.(''MaxOrder'')=1; \n'...
                    SXSTRING{i}(1:end-1) '.(''NumIntSteps'')=2; \n'...
                    SXSTRING{i}(1:end-1) '.(''Length'')=0; \n'...
                    SXSTRING{i}(1:end-1) '.(''Energy'')=' num2str(E0) '; \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomB'')=zeros(1,4); \n'...
                    SXSTRING{i}(1:end-1) '.(''PolynomA'')=zeros(1,4); \n'...
                    ];
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'KICKER','hkicker','vkicker','kicker'}
                
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1) '=atmultipole('...
                    '''' SXSTRING{i}(1:end-1) ''''...
                    ',0,[0 0 0 0],[0 0 0 0],'...
                    '''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''MaxOrder'')=1; \n'...
                    SXSTRING{i}(1:end-1) '.(''Energy'')=' num2str(E0) '; \n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesMADX_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'CONSTAN','constant'}
                % disp('constant')
                
                var=[var SXSTRING{i}(1:end-1) '=' SXSTRING{i+3} '; \n '];
                
                j=j+3;
            case {'sequence'}
                % next element is the length
                seqname=SXSTRING{i}(1:end-1);
                l=SXSTRING{i+j+2};
                j=j+2+1;
                
                elname=SXSTRING{i+j};
                at=SXSTRING{i+j+1};
                lines=[lines seqname '={'];
                sposstring=['\n spos=['];
                %zero=0;
                %driftcounter=1;
                while strcmp(at,'at=') && ~strcmp(elname,'endsequence')
                    
                    elname=SXSTRING{i+j};
                    if ~strcmp(elname,'endsequence')
                        elname=strrep(elname,',','');
                        at=SXSTRING{i+j+1};
                        spos=SXSTRING{i+j+2};
                    
                        lines=[lines upper(elname) ';...\n'];
                        sposstring=[sposstring ',...\n ' num2str(spos)];
            
                        j=j+3;
                    else
                        j=j+1;
                    end
                end
                lines=[lines '};\n'];
                sposstring=[sposstring '];\n'];
                sposstring(10)=[]; % remove extra comma
                % call function that builds sequence from list of elements
                % and total length of sequence
                lines=[lines sposstring '\n %% BUILD LATTICE \n'...
                    seqname '=buildATLattice(' seqname ',spos,' l ');\n'...
                   ... seqname '_red=atreduce(' seqname ');\n'...
                   ... seqname '_FF=FringeSwitch(' seqname ',1);\n'...
                    ];
                
            otherwise
                disp(['Unknown element type: ' ElementType])
        end
        
        
        
    else % variable declaSXSTRING{i}ration??
        if SXSTRING{i}(1)~='!';
           
            % in mad8 declaring a variable before using it is not compulsary.
            % check that all definitions are at the begining.
            if sum(ismember('ABCDFGHILMNOPQRSTUVZWYK',SXSTRING{i+2}))>0 
                % if letters (but E for exponential) then it is a formula
                formulas=[formulas SXSTRING{i} ' = ' SXSTRING{i+1} '; \n '];
            else
                var=[var SXSTRING{i} ' ' SXSTRING{i+1} '; \n '];
            end
            j=2;
        end
    end % if new element
    
    i=i+j;
    j=1;
    
end



%% save close and exit

macroconvertmadXAT=strrep([mst var formulas def lines],';;',';');
fprintf(mafileout,macroconvertmadXAT);

fclose('all');

%% clean workspace
clear macroconvertmadXAT formulas def lines var mst 
clear SXSTRING i j SXCELL elemcount ElementType elname

%% run macro file and save workspace.

[~,seqfilemadX]=fileparts(seqfilemadX);

try % try to run the macro generated
    %%%!!!!! THIS COMAND MAY FAIL! CHECK THE ORDER OF THE DECLARATIONS IN MADX!
    run(filemacroname)

    if nargin<3
        save([seqfilemadX '_AT_LATTICE']);
    else
        save(outfilename);
    end
    
    delete(filemacroname);
    
catch %#ok<CTCH>
    
    fileoutname=[seqfilemadX '_AT_macro.m'];
    
    movefile(filemacroname,fileoutname);
    
    disp(['saved macro in : ' fileoutname]);
    
    error(['could not run the macro file.'...
        ' It is now in your current directory.'...
        ' Please check the order of the definitions']);
end

fclose('all');

return
