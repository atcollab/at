function AT_2_madX(AT_ring,linename)
% function AT_2_madX(AT_ring,linename)
% this functions converts the AT lattice AT_ring in madX format.
% 
% a MADX LINE is generated.
% 
% file ['' linename '_lattice.madx'] is generated contiaining the lattice
% elements definitions and the LINE. no other comands introduced
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

elelat=[elelat '! DEFINITIONS \n\n'];

format='%8.10f';

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
                ' SBEND,L=' num2str(el.('Length'),format) ', '...
                ' ANGLE = ' num2str(el.('BendingAngle'),format) ', '...
                ' K1 = ' num2str(el.('PolynomB')(2),format) ', '...
                ' E1= ' num2str(el.('EntranceAngle'),format) ', '...
                ' E2= ' num2str(el.('ExitAngle'),format) '; '...
                ];
            
            elelat=[elelat di '\n']; %#ok<*AGROW>
        case {'QP','Quadrupole'} % quadrupole
            if el.('MaxOrder')==3
                if el.('PolynomB')(4)==0
                    qp=[' ' el.('FamName') ': '...
                        ' quadrupole,  L = ' num2str(el.('Length'),format)  ', '...
                        'K1 = ' num2str(el.('PolynomB')(2),format) '; '...
                        ];
                else % if octupole in quadrupole, split quadrupole as in mad8 q oc q q oc q
                    nslice=10;
                    qp=[' ' el.('FamName') '_sl: '...
                        ' quadrupole,  L = ' num2str(el.('Length')/nslice/2,format)  ', '...
                        'K1 = ' num2str(el.('PolynomB')(2),format) '; '...
                        '\n'...
                        ];
                    polB=el.('PolynomB');
                    polA=el.('PolynomA');
                    polB(2)=0;
                    polA(2)=0;
                    
                    qp=[qp ' ' el.('FamName')   '_oc: '...
                        ' MULTIPOLE,' ...  LRAD= ' num2str(el.('Length'),format)  ', '...
                        'KNL:= {' num2str(polB.*(factorial(0:length(polB)-1))*el.('Length')/(nslice),[format ', ']) ' 0.0} ,'...
                        'KSL:= {' num2str(polA.*(factorial(0:length(polA)-1))*el.('Length')/(nslice),[format ', ']) ' 0.0};'...
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
                    ' quadrupole,  L = ' num2str(el.('Length'),format)  ', '...
                    'K1 = ' num2str(el.('PolynomB')(2),format) '; '...
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
                ' SEXTUPOLE,  L= ' num2str(el.('Length'),format)  ', '...
                ' K2= ' num2str(el.('PolynomB')(3)*2,format) ' ;'...
                ];
            elelat=[elelat sx '\n'];
        case {'OC','Octupole'} % sextupole
            sx=[el.('FamName')   ': '...
                ' OCTUPOLE,  L= ' num2str(el.('Length'),format)  ', '...
                ' K3=' num2str(el.('PolynomB')(4)*6,format) '*OCT_ON ;'...
                ];
            elelat=[elelat sx '\n'];
        case {'MP','Multipole'} % multipole
            mp=[el.('FamName')   ': '...
                ' MULTIPOLE,  LRAD= ' num2str(el.('Length'),format)  ', '...
                'KNL:= {' num2str(el.('PolynomB').*(factorial(0:length(el.('PolynomB'))-1)),[format ', ']) '0.0} ,'...
                'KSL:= {' num2str(el.('PolynomA').*(factorial(0:length(el.('PolynomA'))-1)),[format ', ']) '0.0};'...
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
            ki=[el.('FamName') ' :KICKER' ' ;'...
                ];
            
            elelat=[elelat ki '\n'];
        case {'DR','Drift'} % drift
            dr=[el.('FamName') ' : DRIFT, L= ' num2str(el.('Length'),format) ' ;'];
            elelat=[elelat dr '\n'];
        case {'CAV','RFCavity'} % drift
            rfc=[el.('FamName') ' : RFCavity, L=' num2str(el.('Length'),format)...
                ',VOLT=RF_ON*'  num2str(el.('Voltage'),format) ' ;'...
                ', freq=' num2str(el.('Frequency'),format) '' ' ;'];
            
            elelat=[elelat rfc '\n'];
            
        otherwise
            warning(['Element: ' el.('FamName') ' was not converted, since it does not match any Class.'])
            mrk=[el.('FamName') ' :MARKER' ' ;'...
                ];
            elelat=[elelat mrk '\n'];
    end
    
end


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


%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return