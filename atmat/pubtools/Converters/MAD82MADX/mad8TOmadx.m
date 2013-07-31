function [seqfilemadX]=mad8TOmadx(seqfilemad8,periodname)
% converts mad8 sequence files to madX
%
%function [seqfileMADX]=mad8TOmadx(seqfilemad8)
%
% This procedure reads a saved sequence in
% mad8 (SAVE,FILE='seqfilemad8';)
% and converts it to madx sequence
% every = goes to :=
% the order of the declarations is the same in the two files.
%
% works also with single mad8 files not containing comands, only
% definitions.
% does not translate call to files since those may change name
% 
% parameters:
%    - seqdilemad8=name of the mad8 lattice file
%    - periodname (optional)= name of the period to use in madx (default is the filename)
% 
% Simone Maria Liuzzo PhD@LNF 25-11-2011
%     update 29-2-2012 : corrected a bug that would not translate correctly
%     markers, kickers and monitor declared only by element name ("BPM: monitor" would not convet properly)

%% open mad8 sequence file and madx sequence output file.
disp('if periodname is not specified I assume the name of the file (up to the .) to be the sequence name that you want to use')
s8=fopen(seqfilemad8,'r');

sxs=[]; % string to input to madx file

% i : segnano tutte le definizioni. ogni : ??? preceduto dal nome
% dell'oggetto e seguito dal tipo di oggetto.
% tra un : e il successivo trovo tutti i parametri di un oggetto e il nome
% dell'oggetto successivo, che sia definizione o linea.


S8CELL=textscan(s8,'%s');
S8STRING=S8CELL{1}; % still a cell array but every space or new line now is stored as a new cell.


%save('1');

%return
%% convert to madx
j=1; % used in line block counter
sxs=['! MADX converted sequence from ' seqfilemad8 ' '];
corrcount=1;
i=12; % skip header
S8STRING{12}
while i<length(S8STRING)-2
    
    if S8STRING{i}(end)==':' % new element or line
        sxs=[sxs ';\n ' S8STRING{i}]; %#ok<*AGROW> % end line and go to newline add name
       
        S8STRING{i+j}=strrep(S8STRING{i+j},'MARKER','MARKER,');
        S8STRING{i+j}=strrep(S8STRING{i+j},'MARKER,,','MARKER,');
        S8STRING{i+j}=strrep(S8STRING{i+j},'KICKER','KICKER,');
        S8STRING{i+j}=strrep(S8STRING{i+j},'KICKER,,','KICKER,');
        S8STRING{i+j}=strrep(S8STRING{i+j},'MONITOR','MONITOR,');
        S8STRING{i+j}=strrep(S8STRING{i+j},'MONITOR,,','MONITOR,');

        ElementType=S8STRING{i+j}(1:end-1);
        %  sxs=[sxs ' ' S8STRING{i+1}];
        %j=2;
        nwel=S8STRING{i+j}(end);
        switch ElementType
            case 'QUADRUPOLE'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
                % build in monitor, and corrector
                    elname=S8STRING{i}(1:end-1);
                     if strcmp(elname(1:2),'QD') 
                     else
                      Q_CM=[';\n'...
                        'CH_' num2str(corrcount) ': Hkicker, kick:=0,L=0;\n'...
                        'CV_' num2str(corrcount) ': Vkicker, kick:=0,L=0;\n'...
                        'MH_' num2str(corrcount) ': Hmonitor;\n'...
                        'MV_' num2str(corrcount) ': Vmonitor;\n'...
                        'h_' num2str(corrcount) ': QUADRUPOLE, K1=' elname '->K1, L=(' elname '->L)/2;\n'...
                        elname ': LINE=(h_' num2str(corrcount) ',CH_' num2str(corrcount) ',CV_' num2str(corrcount) ',MH_' num2str(corrcount) ', MV_' num2str(corrcount) ',h_' num2str(corrcount) ');\n'...
                        ];
                    corrcount=corrcount+1;

                   % sxs=[sxs ' ' Q_CM];
                    end     
            case 'SEXTUPOLE'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
                
                    elname=S8STRING{i}(1:end-1);
                    if strcmp(elname(1:2),'SD') 
                    S_CM=[';\n'...
                        'CH_' num2str(corrcount) ': Hkicker, kick:=0,L=0;\n'...
                        'CV_' num2str(corrcount) ': Vkicker, kick:=0,L=0;\n'...
                        'ksq' num2str(corrcount) ':=0;\n'...
                        'SQ_' num2str(corrcount) ': MULTIPOLE, ksl:={0,ksq' num2str(corrcount) '},L=0;\n'...
                        'MH_' num2str(corrcount) ': Hmonitor;\n'...
                        'MV_' num2str(corrcount) ': Vmonitor;\n'...
                        'h_' num2str(corrcount) ': SEXTUPOLE, K2:=' elname '->K2, L:=(' elname '->L)/2;\n'...
                        elname ': LINE=(h_' num2str(corrcount) ',CH_' num2str(corrcount) ',CV_' num2str(corrcount) ',SQ_' num2str(corrcount) ',MH_' num2str(corrcount) ', MV_' num2str(corrcount) ',h_' num2str(corrcount) ');\n'...
                        ];
                    
                   % sxs=[sxs ' ' S_CM];
                  corrcount=corrcount+1;                       
                    end     
            case 'RBEND'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'SBEND'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'DRIFT'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'RFCAVITY'
                while nwel~=':'
                    a=sum(strfind(S8STRING{i+j},'REMOVED'));
                    if a==0
                        
                        sxs=[sxs ' ' S8STRING{i+j}];
                    else
                        disp(S8STRING{i+j})
                    end
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                    
                end
            case 'MULTIPOLE'  %%%%% CHANGE TO GOOD MULTIPOLE!
                KNL=['KNL:={'];
                while nwel~=':'
                    if S8STRING{i+j}(1)~='K'
                        sxs=[sxs ' ' S8STRING{i+j}];
                    else
                        for knl=1:str2double(S8STRING{i+j}(2))
                            KNL=[KNL '0,'];
                        end
                        KNL=[KNL S8STRING{i+j}(5:end) '*OCT_ON}'];
                        sxs=[sxs ' ' KNL ];
                    end
                    %        sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'OCTUPOLE'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'SOLENOID'
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end

            case 'MARKER'  %%%%% CHANGE TO GOOD MULTIPOLE!
                while nwel~=':'
                    
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'MONITOR'  %%%%% CHANGE TO GOOD MULTIPOLE!
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'KICKER'  %%%%% CHANGE TO GOOD MULTIPOLE!
                while nwel~=':'
                    sxs=[sxs ' ' S8STRING{i+j}];
                    j=j+1;
                    nwel=S8STRING{i+j}(end);
                end
            case 'CONSTAN'  %%%%% CHANGE TO GOOD MULTIPOLE!
                sxs=[S8STRING{i}(1:end-1) '=' S8STRING{i+j+2} ';\n' sxs(1:end-length(S8STRING{i+j}))];
                j=j+3;
            otherwise
                %%% unrecognized LINE comand
                LineMAD=[];
                if strcmp(S8STRING{i+j}(1:5),'LINE=')
                    nwel2=S8STRING{i+j+1}(end);
                    while nwel~=':' && nwel2~='=' %second case is the value definition
                        % clean from &
                        if S8STRING{i+j}(end)=='&';
                            S8STRING{i+j}=S8STRING{i+j}(1:end-1);
                        end
                        
                        % apend to madx file
                        LineMAD=[LineMAD ' ' S8STRING{i+j}];
                        j=j+1;
                        nwel=S8STRING{i+j}(end);
                        nwel2=S8STRING{i+j+1}(end);
                    end
                    % transform ,n*QDF, in ,n*(QDF),
                    multiplication=strfind(LineMAD,'*');
                    k=0;
                    for im=multiplication
                        ind=im+k;
                        comma=strfind(LineMAD(ind:end),',');
                        
                        if length(comma)>=1 && ~isempty(LineMAD(ind:end))
                            comma=comma(1); % get first comma after multiplication
                            LineMAD=[LineMAD(1:ind) '(' LineMAD(ind+1:ind+comma-2) ')' LineMAD(ind+comma-1:end)];
                        else
                            LineMAD=[LineMAD(1:ind) '(' LineMAD(ind+1:end) ')'];
                        end
                        k=k+2; % length is encreased by 2 units at every loop
                    end
                    sxs=[sxs ' ' LineMAD ';\n'];
                    
                    if  nwel=='='
                        j=j-1;
                        %                    S8STRING{i}
                        %                    j
                        %                    S8STRING{i+j}
                    end
                    
                else
                    disp(['Unknown element type: ' ElementType])
                end
                
                
        end
        
        
        
    else % variable declaration??
        
        % copy 3 blocks on top of all variables
        sxs=[S8STRING{i} S8STRING{i+1} S8STRING{i+2} ';\n ' sxs];
        j=3;
    end % if new element
    
    i=i+j;
    j=1;
    
end

disp(corrcount)

%% run madx sequence file to check tunes, length betas and dispersion.

%% out madX sequence file
[~,seqfilemad8]=fileparts(seqfilemad8);
seqfilemadX=[seqfilemad8 'X'];
sX=fopen(seqfilemadX,'w+');
sXt=fopen(['test_' seqfilemadX ],'w+');


if (nargin == 1) % only 1 argument provided
    sqname=strrep(seqfilemadX,'.mad8X','');
    sqname=strrep(sqname,'.madX','');
    sqname=strrep(seqfilemadX,'.MAD8X','');
    sqname=strrep(sqname,'.MADX','');
    sxs=strrep(sxs,'RING_FF',sqname);
elseif (nargin == 2) % provide also name of period to use in test
    
    sqname=periodname;
    
    sxs=strrep(sxs,'RING_FF',sqname);
end

testtwiss=[...
'setplot,ascale=2,lscale=2,sscale=1.5,rscale=2,xsize=28,ysize=20, \n'...
 '       post=2,font=-4,lwidth=10;\n'...
    '\n beam;'...
    '\n use, period=' sqname ';'...
    '\n twiss,chrom,file=test.twiss;'...
    '\n plot, haxis=s, vaxis=rbetx,rbety, interpolate, style=1,colour=100,file=' sqname ';'...
...     '\n plot, haxis=s, vaxis=Wx,Wy, interpolate,style=1,colour=100;'...
...     '\n plot, haxis=s, vaxis=dx,dy, interpolate,style=1,colour=100;'...
...     '\n plot, haxis=s, vaxis=x,y, interpolate,style=1,colour=100;'...
...     '\n plot, range=#s/SF1L[2], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n plot, range=SF1L[4]/SF1L[8], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n plot, range=SF1L[8]/SF1L[12], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n plot, range=SF1L[12]/SF1[2], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n plot, range=SF1[2]/SF1[8], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy,  interpolate,style=1,colour=100;'...
...     '\n plot, range=SF1[8]/SF1[12], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy,  interpolate,style=1,colour=100;'...
...     '\n plot, range=SF1[12]/#e, haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n plot, range=OC1[1]/CRABL[2], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;'...
...     '\n  plot, range=CRABL[2]/IP[1], haxis=s, vaxis1=betx,bety, vaxis2=dx,dy, interpolate, style=1,colour=100;\n'...
     ];

trackingLines=[...
    '\n select, flag=twiss,class=Hkicker,column=name,s; twiss,file=tt;\n\n'...
    '\n SURVEY, x0=0, y0=0, z0=0, theta0=0, phi0=0, psi0=0,file=surv, table=survey, sequence=sqname;\n'...
    'plot, table=survey,haxis=z,vaxis=x;\n'...,hmin=-250,hmax=250,vmin=-450,vmax=50;\n'...
    'save, sequence=' sqname ', file=' seqfilemadX(~ismember(seqfilemadX,'.')) '.seq; ! savesequence \n '...
    ...'use,period=ring_noff;\nsave, sequence=RING_noff, file=' seqfilemadX(~ismember(seqfilemadX,'.')) 'noff.seq,bare; ! savesequence \n '...
    ... '    savebeta,label=HER,place=#e;\n'...
    ];

rootdrawlayout=[
'\n use, period=' sqname ';\n'...
'ptc_create_universe;\n'...
'ptc_create_layout, model=1, method=6, nst=100, exact=false, closed_layout=false;\n'...
'ptc_setswitch, debuglevel=2, maxacceleration=true, exact_mis=true, time=true, totalpath=true, fringe=true;\n'...
'ptc_printframes, file="draw' sqname '.C", format=rootmacro;\n'...
'ptc_end;\n'...
];


sxs=strrep(sxs,'&',''); % remove &
sxs=strrep(sxs,';\n;',';\n\n'); % remove useless ;
sxs=strrep(sxs,'K2=','K2=SXT_ON*');  % add sextupole on switch
sxs=strrep(sxs,'=',':='); % add :=
sxs=strrep(sxs,'::=',':=');
sxs=strrep(sxs,']',''); % convert [] to ->
sxs=strrep(sxs,'[','->');
sxs=strrep(sxs,'LINE:=','LINE=');  % restore Line definitions to
sxs=strrep(sxs,', ',',');


%% save close and exit
fprintf(sX,[sxs]);
fprintf(sXt,['call,file=' seqfilemadX ';\n' testtwiss trackingLines rootdrawlayout]);

system(['madxp < ' ['test_' seqfilemadX ] ])

fclose('all');
clear all
close all

return