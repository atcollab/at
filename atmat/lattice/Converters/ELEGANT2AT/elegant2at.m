function elegant2at(elegantlattice,E0,outfilename)
%function elegant2at(elegantlattice,E0,outfilename)
% tansform elegant %s.new file (save_lattice with output_seq=0) file into AT lattice structure.
%
% This procedure reads a saved elegant lattice and converts it to an AT lattice
%
% Elegant command to save the lattice sequence : 
% 
%  _______ filename.ele_________
%  &save_lattice
%   filename = %s.new,
%   output_seq=0,
% &end
%  ___________________________
%
%  filename.new will contain the list of all the elements that make up the
%  lattice in the correct format in a single file
% 
% The routine outputs a Matlab macro with all the AT defitions and variables as
% in the elegant file  
%
% Works also with single elegant file not containing commands, only
% definitions.
%
% parameters:
%    - elegantlattice = name of the elegant lattice file
%    - E0  = design energy
%    - outfilename (default: elegantlattice_AT_LATTICE.mat)
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
% created 7-sep-2012 by S.M.Liuzzo @ ESRF as atfrommadx
% modified 1-12-2015 by G. Campogiani @LNF as elegant2at
                  

%% initial warnings
disp('PRELIMINARY INFORMATION:')
disp('');
disp(['1) THE elegant FILE MUST BE OUTPUT OF save_lattice', ...
   '&save_lattice filname = %s.new, output_seq=0, &end'...
    ' THIS GARANTES AN APPROPRIATE FILE FORMAT']);
disp('');
disp(['2) THE elegant PROGRAM allows to use variable not previously defined.'...
    ' This is not known to AT!.'...
    ' If the program fails but generates a ..._macro.m file, '...
    ' please edit this file reordering the declarations and run it.']);

%% open madX sequence file

% the : detect definitions. every : is preceded by a name
% of an element.
% between two : all the parameter of an elemetn or a line ar found.
sX=fopen(elegantlattice,'r');

%% get madX file in a cell array with a new elment at every space comma or new line.

SXCELL=textscan(sX,'%s','delimiter','\n');
% SXCELL=strrep(cellstr(SXCELL),'\n','&');

SXSTRING=SXCELL{1}; % still a cell array but every space or new line now is stored as a new cell.

% scroll and reshape to divide atributes
tmp={};

for i=1:length(SXSTRING)
    if SXSTRING{i}(end)=='&'
        SXSTRING{i}=SXSTRING{i}(1:end-1);
    end;

    SXSTRING{i}=strrep(SXSTRING{i},' ',''); % no spaces



    
    c=[1 sort([strfind(SXSTRING{i},',') strfind(SXSTRING{i},':')...
                     strfind(SXSTRING{i},'=') strfind(SXSTRING{i},'(')...
                     ]) length(SXSTRING{i})];
    
    for jc=1:length(c)
        if jc==1
            tmp=[tmp SXSTRING{i}(c(jc):c(jc+1))];
        else
            if jc==length(c)
                tmp=[tmp SXSTRING{i}(c(jc))];
            else
                tmp=[tmp SXSTRING{i}(c(jc)+1:c(jc+1))];
            end
        end
    end

    
end
SXSTRING=tmp;



%% open .m file to output matlab translated code

filemacroname=[tempname '.m'];

mafileout=fopen(filemacroname,'w+');
mst=['%% this is a macro that converts to AT the madX lattice: '...
                     elegantlattice '\n%%\n%% Created: ' datestr(now)...
                     '\n%%\n%%\n\nglobal GLOBVAL;\nGLOBVAL.E0=' num2str(E0)...
                     ';\n\n\n'];
def=['\n\n%%%% DEFINITIONS \n\n']; %#ok<*NBRAK>
lines=['\n\n%%%% LINES \n\n'];

%% convert to a matlab macro
j=1; % used in line block counter (element attributes counter)

elemcount=0;
i=1; 

while i<length(SXSTRING)-2
    
    if SXSTRING{i}(end)==':' % new element or line
        def=[ def ...
            ]; % end line and go to newline add name
        SXSTRING{i}(1:end-1)=upper(SXSTRING{i}(1:end-1));

%          not very clear what's the purpose of the following lines        
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'MARKER','MARKER,');
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'MARKER,,','MARKER,');
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'KICKER','KICKER,');
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'KICKER,,','KICKER,');
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'MONITOR','MONITOR,');
%         SXSTRING{i+j}=strrep(SXSTRING{i+j},'MONITOR,,','MONITOR,');
        
        ElementType=SXSTRING{i+j}(1:end-1);

        nwel=SXSTRING{i+j}(end);
      
        switch ElementType
            case {'QUAD','quadrupole','QUADRUPOLE', 'QUADRUPO'}
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atquadrupole('''...
                    SXSTRING{i}(1:end-1)...
                    ''',0,0,''StrMPoleSymplectic4Pass'');\n'...
                    SXSTRING{i}(1:end-1) '.(''NumIntSteps'')=10; \n'...
                    ];
                
                
                elemcount=elemcount+1;
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
            case {'SEXT','sextupole','SEXTUPOLE'}
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atsextupole('''...
                    SXSTRING{i}(1:end-1)...
                    ''',0,0,''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                    
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
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),SXSTRING{i+j},SXSTRING{i+j+1});
                    
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
              
            case {'SBEN','sbend','SBEND'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atsbend('...
                    '''' SXSTRING{i}(1:end-1)...
                    ''',0,0,0,''BndMPoleSymplectic4Pass'');\n'...
                    SXSTRING{i}(1:end-1) '.(''NumIntSteps'')=10; \n'...
                    ];
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
           def=[ def SXSTRING{i}(1:end-1) '.(''MaxOrder'')=length(' SXSTRING{i}(1:end-1) '.(''PolynomB''))-1; \n'];
              
            case {'DRIF','DRIFT','drift'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atdrift('...
                    '''' SXSTRING{i}(1:end-1)...
                    ''',0,''DriftPass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'RFC','RFCA','RFCAVITY','rfcavity'}
                def=[def '\n'];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1)...
                    '=atrfcavity(''' SXSTRING{i}(1:end-1)...
                    ''',0,0,0,0,' num2str(E0) ...
                    ',''DriftPass'');\n'...
                    ];
                
                elemcount=elemcount+1;
                
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'MULT','MULTIPOLE','multipole'}
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
                    
                    if ~isempty(multipoles)
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},multipoles);
                
                    else
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                    SXSTRING{i+j},SXSTRING{i+j+1});
                    end
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
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
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
                    '''SolenoidLinearPass'');\n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
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
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
            case {'HMO','VMO','HMON','VMON'}
                def=[def '\n'];
%                      SXSTRING{i}(1:end-1) '.(''FamName'')=''' SXSTRING{i}(1:end-1) ''';\n' ...
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''FamName'')=''BPM'';\n' ...
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
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
            case {'MAR','MARK'}
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''FamName'')=''' SXSTRING{i}(1:end-1) ''';\n' ...
                    SXSTRING{i}(1:end-1) '.(''Class'')=''Monitor''; \n'...
                    SXSTRING{i}(1:end-1) '.(''Length'')=0; \n'...
                    SXSTRING{i}(1:end-1) '.(''PassMethod'')=''IdentityPass''; \n'...
                    ];
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end

            case {'HKIC','VKIC','KIC','hkic','vkic','kic'}
                
                def=[def '\n'];
                def=[ def ...
                    SXSTRING{i}(1:end-1) '= atmultipole('...
                    '''' SXSTRING{i}(1:end-2) ''''...
                    ',0,[0 0 0 0],[0 0 0 0],'...
                    '''StrMPoleSymplectic4Pass'');\n'...
                    ];
                
                def=[ def ...
                    SXSTRING{i}(1:end-1) '.(''MaxOrder'')=1; \n'...
                    SXSTRING{i}(1:end-1) '.(''Energy'')=' num2str(E0) '; \n'...
                    SXSTRING{i}(1:end-1) '.(''NumIntSteps'')=10; \n'...
                    ];
                
                
                elemcount=elemcount+1;
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
            case {'MAXAMP','maxamp'}
                def = [def SXSTRING{i}(1:end-1) '= ataperture('...
                    '''' SXSTRING{i}(1:end-1) ''''...
                    ',[0 0 0 0],AperturePass'');\n'...
                    ];
                
                while nwel~=':' % loops atributes of this element definition
                    def=ParseAtributesELEGANT_2_AT(def,SXSTRING{i}(1:end-1),...
                        SXSTRING{i+j},SXSTRING{i+j+1});
                    
                    j=j+1; %go to new atribute
                    nwel=SXSTRING{i+j}(end);
                end
                
            case {'LINE','line'}
                % next element is the length
                seqname=SXSTRING{i}(1:end-1);
                j=j+1;
                
                elname=SXSTRING{i+j};
                at='at =';
                lines=[lines seqname '= {'];
                
                while strcmp(at,'at =') && ~strcmp(elname,')')
                elname= SXSTRING{i+j};

                    if ~strcmp(elname,')')
                        elname=strrep(elname,' ','');
                        elname=strrep(elname,'  ','');
                        elname=strrep(elname,',','');
                        elname=strrep(elname,'(','');
                        elname=strrep(elname,')','');
                
                        lines=[lines upper(elname) ' '];
            
                        j=j+1;
                        
                    end                   
                end
                lines=[lines '};\n'];

%                 sposstring=[sposstring '];\n'];
%                 sposstring(10)=[]; % remove extra comma
                % call function that builds sequence from list of elements
                % and total length of sequence
                
            otherwise
                disp(['Unknown element type: ' ElementType])
            end
    end % if new element
    
    i=i+j;
    j=1;
    
end
lines = [lines '\n %% BUILD LATTICE \n'...
    '\nglobal THERING;\nTHERING = transpose(' seqname ');\n'];


%% save close and exit

macroconvertmadXAT=strrep([mst def lines],';;',';');
fprintf(mafileout,macroconvertmadXAT);

fclose('all');

%% clean workspace
clear macroconvertmadXAT def lines mst 
clear SXSTRING i j SXCELL elemcount ElementType elname

%% run macro file and save workspace.

[~,elegantlattice] = fileparts(elegantlattice);

try % try to run the macro generated
    %%%!!!!! THIS COMAND MAY FAIL! CHECK THE ORDER OF THE DECLARATIONS IN MADX!
    run(filemacroname)
        
    if nargin<3
        fileoutname = [elegantlattice '_AT.m'];
    end
    
    name=strrep(fileoutname,'.m','');
    varoutname=sprintf('%s.mat',name);  

    save(varoutname,'THERING');
    
    fclose('all');

    movefile(filemacroname,fileoutname);
%   delete(filemacroname);
    

catch %#ok<CTCH>
    
    fileoutname=[elegantlattice '_AT_macro.m'];
    movefile(filemacroname,fileoutname);

    disp(['saved macro in : ' fileoutname]);
    
    error(['Could not run the macro file.'...
        ' It is now in your current directory.'...
        ' Please check the order of the definitions']);
end


return
