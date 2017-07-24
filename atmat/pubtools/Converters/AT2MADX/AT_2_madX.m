function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)
% function [elelat,defs,lines]=AT_2_madX(AT_ring,linename)
% this functions converts the AT lattice AT_ring in madX format.
% 
% a MADX LINE is generated.
% 
% file ['' linename '_lattice.madx'] is generated contiaining the lattice
% (elelat)  elements definitions (defs) and the LINE (lines). 
% no other comands are introduced
%
% to test in MADX run (replace linename with apropriate name)
% madx < linename_lattice.madx 
% 
% to test with twiss and plot add following lines at the end of file:
% beam;
% use, period=linename;
% twiss;
% plot, haxis=s, vaxis1=betx,bety, vaxis2=dx;
% 

outfile=['' linename '_lattice.madx'];

%outfile='madXelemdef.elem';
elelat=['!!\n!!  madX lattice elements: ' linename '\n!!  Created: ' datestr(now) '\n!!\n!!\n\n'];

%% get family names for definitions
[families,ind_first_oc_ring]=...
    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');
formatvalues='%10.12f';

elelat=[elelat '! DEFINITIONS \n\n BEAM, PARTICLE=ELECTRON, ENERGY=' num2str(atenergy(AT_ring)*1e-9,formatvalues) ';\n OCT_ON:=1;\n SXT_ON:=1;\n RF_ON:=0;'];


%% loop families for definitions
for i=1:length(families)
    el= AT_ring{ind_first_oc_ring(i)};
    
    if isfield(el,'BetaCode')
        type=el.BetaCode;
    elseif isfield(el,'Class')
        type=el.Class;
    else
        type='Marker';
    end
    
    switch type
        case {'DI','Dipole','Bend'} % dipole
            
            di=[el.('FamName')   ': '...
                ' SBEND,L=' num2str(el.('Length'),formatvalues) ', '...
                ' ANGLE = ' num2str(el.('BendingAngle'),formatvalues) ', '...
                ' K0 = ' num2str(el.('PolynomB')(1),formatvalues) ', '...
                ' K1 = ' num2str(el.('PolynomB')(2),formatvalues) ', '...
                ' E1= ' num2str(el.('EntranceAngle'),formatvalues) ', '...
                ' E2= ' num2str(el.('ExitAngle'),formatvalues) ', '...
                ];
            
            if isfield(el,'FringeInt1')
            di=[di ' FINT= ' num2str(el.('FringeInt1'),formatvalues) ', '...
                ' FINTX= ' num2str(el.('FringeInt2'),formatvalues) ', '...
                ' hgap= ' num2str(el.('FullGap'),formatvalues) '/2, '...
                ];
            end
            
            di=[di ';'];
            
            elelat=[elelat di '\n']; %#ok<*AGROW>
        case {'QP','Quadrupole'} % quadrupole
            if el.('MaxOrder')==3
                if el.('PolynomB')(4)==0
                    qp=[' ' el.('FamName') ': '...
                        ' quadrupole,  L = ' num2str(el.('Length'),formatvalues)  ', '...
                        'K1 = ' num2str(el.('PolynomB')(2),formatvalues) '; '...
                        ];
                else % if octupole in quadrupole, split quadrupole as in mad8 q oc q q oc q
                    nslice=10;
                    qp=[' ' el.('FamName') '_sl: '...
                        ' quadrupole,  L = ' num2str(el.('Length')/nslice/2,formatvalues)  ', '...
                        'K1 = ' num2str(el.('PolynomB')(2),formatvalues) '; '...
                        '\n'...
                        ];
                    polB=el.('PolynomB');
                    polA=el.('PolynomA');
                    polB(2)=0;
                    polA(2)=0;
                    
                    qp=[qp ' ' el.('FamName')   '_oc: '...
                        ' MULTIPOLE,' ...  LRAD= ' num2str(el.('Length'),format)  ', '...
                        'KNL:= {' num2str(polB.*(factorial(0:length(polB)-1))*el.('Length')/(nslice),[formatvalues ', ']) ' 0.0} ,'...
                        'KSL:= {' num2str(polA.*(factorial(0:length(polA)-1))*el.('Length')/(nslice),[formatvalues ', ']) ' 0.0};'...
                        '\n'...
                        ];
                    
                    qp=[qp el.('FamName') ': line=('...
                        repmat([el.('FamName') '_sl,'...
                        el.('FamName') '_oc,'...
                        el.('FamName') '_sl,'...
                        ],1,nslice)...
                        '\n'...
                        ];
                    qp(end-2:end)=[];
                    qp=[qp ');\n'];
                end
            else
                qp=[' ' el.('FamName') ': '...
                    ' quadrupole,  L = ' num2str(el.('Length'),formatvalues)  ', '...
                    'K1 = ' num2str(el.('PolynomB')(2),formatvalues) '; '...
                    ];
                
            end
            
            
            elelat=[elelat qp '\n\n'];
            %             qp=[el.('FamName')   ': '...
            %                 ' QUADRUPOLE,  L= ' num2str(el.('Length'),format)  ','...
            %                 ' K1= ' num2str(el.('PolynomB')(2),format) ';'...
            %                 ];
            %
            %             elelat=[elelat qp '\n'];
        case {'SX','Sextupole'} % sextupole
            sx=[el.('FamName')   ': '...
                ' SEXTUPOLE,  L= ' num2str(el.('Length'),formatvalues)  ', '...
                ' K2= ' num2str(el.('PolynomB')(3)*2,formatvalues) '*SXT_ON ;'...
                ];
            elelat=[elelat sx '\n'];
        case {'OC','Octupole'} % sextupole
            sx=[el.('FamName')   ': '...
                ' OCTUPOLE,  L= ' num2str(el.('Length'),formatvalues)  ', '...
                ' K3=' num2str(el.('PolynomB')(4)*6,formatvalues) '*OCT_ON ;'...
                ];
            elelat=[elelat sx '\n'];
        case {'MP','Multipole'} % multipole
            el.('PolynomB')(isnan(el.('PolynomB')))=0.0;
            if el.Length==0
            mp=[el.('FamName')   ': '...
                ' MULTIPOLE,  LRAD= ' num2str(el.('Length'),formatvalues)  ', '...
                ' L= ' num2str(el.('Length'),formatvalues)  ', '...
                'KNL:= {' num2str(el.('PolynomB').*(factorial(0:length(el.('PolynomB'))-1)),[formatvalues ', ']) '0.0} ,'...
                'KSL:= {' num2str(el.('PolynomA').*(factorial(0:length(el.('PolynomA'))-1)),[formatvalues ', ']) '0.0};'...
                ];
            else
            mp=[el.('FamName')   ': '...
                ' Drift,  L= ' num2str(el.('Length'),formatvalues)  '; '...
                ];
            warning('Multipoles have zero length in MADX, CONVERTED TO DRIFT. please model in AT as Drift-lense-Drift.')
            end
            elelat=[elelat mp '\n'];
        case {'ThinMultipole'} % multipole
            formatm='%10.5e';
            mp=[el.('FamName')   ': '...
                ' MULTIPOLE,  LRAD= ' num2str(el.('Length'),formatvalues)  ', '...
                'KNL:= {' num2str(el.('PolynomB').*(factorial(0:length(el.('PolynomB'))-1)),[formatm ', ']) '0.0} ,'...
                'KSL:= {' num2str(el.('PolynomA').*(factorial(0:length(el.('PolynomA'))-1)),[formatm ', ']) '0.0};'...
                ];
            elelat=[elelat mp '\n'];
        case {'PU','Monitor'} % bpm
            pu=[el.('FamName') ' :MONITOR' ' ;'...
                ];
            elelat=[elelat pu '\n'];
        case {'MRK','Marker'}% marker
            mrk=[el.('FamName') ' :MARKER' ' ;'...
                ];
            elelat=[elelat mrk '\n'];
        case {'KI','Corrector'} % kicker
            ki=[el.('FamName') ' :KICKER, L= ' num2str(el.('Length'),formatvalues)  ' '...
                ', HKick:=' num2str(el.PolynomB(1))...
                ', VKick:=' num2str(el.PolynomA(1)) ' ;'...
                ];
            
            elelat=[elelat ki '\n'];
        case {'SKW','SkewQuadrupole'} % kicker
            skw=[el.('FamName') ' :Multipole, KSL:={0,' num2str(el.PolynomA(2))...
                 '} ;'...
                ];
            
            elelat=[elelat skw '\n'];
        case {'DR','Drift'} % drift
            dr=[el.('FamName') ' : DRIFT, L= ' num2str(el.('Length'),formatvalues) ' ;'];
            elelat=[elelat dr '\n'];
        case {'CAV','RFCavity'} % drift
            rfc=[el.('FamName') ' : RFCavity, L=' num2str(el.('Length'),formatvalues)...
                ',VOLT=RF_ON*'  num2str(el.('Voltage')/1e6,formatvalues) ' '...% MV
                ...', NO_CAVITY_TOTALPATH = false'...
                ', freq=' num2str(el.('Frequency')*1e-6,formatvalues) ''... MHz
                ' ;'];
            
            elelat=[elelat rfc '\n'];
            
        otherwise
            warning(['Element: ' el.('FamName') ' was not converted (marker), since it does not match any Class.'])
            elelat=[elelat el.('FamName') ' : marker; \n'];
           
            mrk=[el.('FamName') ' :MARKER' ' ;'...
                ];
            elelat=[elelat mrk '\n'];
    end
    
end
defs=elelat;

elelat=[];

elelat=[elelat '! LINE \n\n'];

elelat=[elelat linename ' : LINE = (' ' '];

%% define lattice line
% loop all elements
for i=1:length(AT_ring)
    if i~=1
        elelat=[elelat ',' AT_ring{i}.('FamName') '\n '];
    else
        elelat=[elelat ' ' AT_ring{i}.('FamName') '\n '];
    end
end

elelat=[elelat ') ;'];
lines=elelat;

elelat=[defs lines];

%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return