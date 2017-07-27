function [elelat,def,lines]=AT_2_mad8(AT_ring,linename)
% function [elelat,def,lines]=AT_2_mad8(AT_ring,linename)
% this functions converts the AT lattice AT_ring in mad8 format.
% 
% a MAD8 LINE is generated.
% 
% file ['' linename '_lattice.mad8'] is generated contiaining the lattice
% (elelat) elements definitions (defs) and the LINE (lines). no other comands introduced
% 
% 

outfile=['' linename '_lattice.mad8'];

elelat=['!!\n!!  mad8 lattice: ' linename '\n!!  Created: ' datestr(now) '\n!!\n!!\n\n'...
    ...'assign,print="' linename 'sliced.print"\n'...
...'assign,echo="' linename 'sliced.echo"\n'...
...'option, -echo\n'...
];

%% get family names for definitions
[families,ind_first_oc_ring]=...
    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');

elelat=[elelat '! DEFINITIONS \n\n'];

format='%8.8f';

%% loop families for definitions
for i=1:length(families)
   el= AT_ring{ind_first_oc_ring(i)};
   if isfield(el,'Class')
    switch el.('Class')
        case 'Bend' % dipole
            %Ang=el.('BendingAngle');
            %disp('all dipoles are considered rbend')
            di=[el.('FamName')   ': '...*sin(Ang/2)/Ang/2
                ' SBEND,L=' num2str(el.('Length'),format) ', '...&\n
                ' ANGLE = ' num2str(el.('BendingAngle'),format) ', &\n'...
                ' K1 = ' num2str(el.('PolynomB')(2),format) ' '...&\n
                ' E1= ' num2str(el.('EntranceAngle'),format) ', '...&\n
                ' E2= ' num2str(el.('ExitAngle'),format) ' ;'...\n\n
               ' '];
           
            elelat=[elelat di '\n']; %#ok<*AGROW>
        case 'Quadrupole' % quadrupole
            
            if strfind(el.PassMethod,'Fringe') 
            qp=[el.('FamName')   ': '...&\n
                ' QUADRUPOLE,  L= ' num2str(el.('Length'),format)  ', '...&\n
                ' K1= ' num2str(el.('PolynomB')(2),format) ',FRGF=.TRUE.;'...\n
                ];
            else
            qp=[el.('FamName')   ': '...&\n
                ' QUADRUPOLE,  L= ' num2str(el.('Length'),format)  ', '...&\n
                ' K1= ' num2str(el.('PolynomB')(2),format) ';'...\n
                ];
            end
            
            elelat=[elelat qp '\n'];
        case 'Sextupole' % sextupole
            sx=[el.('FamName')   ': '...&\n
                ' Sextupole,  L= ' num2str(el.('Length'),format)  ','... &\n
                ' K2= ' num2str(el.('PolynomB')(3)*2,format) ' ;'...\n
                ];
            elelat=[elelat sx '\n'];
        case 'Multipole' % bpm
            % el.('PolynomB')
            ord=find(el.('PolynomB'));
            if isempty(ord)
                ord=length(el.('PolynomB'));
            end
                mp=[el.('FamName')   ': '...&\n
                    ' Multipole, '.... L= ' num2str(el.('Length'),format)  ', &\n'...
                    ...'ORDER= ' num2str(ord(1)-1) ', &\n' ... 
                   'K' num2str(ord(1)-1) 'L= ' num2str(el.('PolynomB')(ord(1))*factorial(ord(1)-1),format) '; \n'...
 %                   'KNL= ' num2str(factorial(ord(1)-1)*el.('PolynomB')(ord(1)),format) ', \n'...
                %' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL+1),'%d') ' \n'...
                ];
            elelat=[elelat mp '\n'];
        case {'ThinMultipole'} % bpm
            formatm='%10.5e';
            % el.('PolynomB')
            ord=find(el.('PolynomB'));
            if isempty(ord)
                ord=length(el.('PolynomB'));
            end
                mp=[el.('FamName')   '_n: '...&\n
                    ' Multipole, '.... L= ' num2str(el.('Length'),format)  ', &\n'...
                    ...'ORDER= ' num2str(ord(1)-1) ', &\n' ... 
                   ];
               
               if ord>9, ord=9; end; % mad8 limit
               
               for imult=1:ord
                   mp=[mp 'K' num2str(imult-1) 'L= ' num2str(el.('PolynomB')(imult)*el.('Length')*factorial(imult-1),format) ', '...
                   ];
               end
               
               mp=[mp(1:end-2) '; \n'];
               
               % skew
               mpa=[el.('FamName')   '_s: '...&\n
                    ' Multipole, '.... L= ' num2str(el.('Length'),format)  ', &\n'...
                    ...'ORDER= ' num2str(ord(1)-1) ', &\n' ... 
                   ];
               
               if ord>9, ord=9; end; % mad8 limit
               
               for imult=1:ord
                   mpa=[mpa 'K' num2str(imult-1) 'L= '...
                       num2str(el.('PolynomB')(imult)*el.('Length')*factorial(imult-1),format) ', '...
                       'T' num2str(imult-1) '= 3.142128/' num2str(2*imult+2) ', '...
                   ];
               end
               
               mpa=[mpa(1:end-2) '; \n'];

               % multipole skew and normal to single multipole element
               m=[el.('FamName')   ': LINE(' el.('FamName') '_n,' el.('FamName') '_s' ');\n'];
               
            elelat=[elelat mp mpa m '\n'];
        case 'Octupole' % bpm
            % el.('PolynomB')
            ord=find(el.('PolynomB'));
            if isempty(ord)
                ord=length(el.('PolynomB'));
            end
                mp=[el.('FamName')   ': '...&\n
                    ' Octupole, L= ' num2str(el.('Length'),format)  ', '...&\n
                    ...'ORDER= ' num2str(ord(1)-1) ', &\n' ... 
                   'K' num2str(ord(1)-1) '= ' num2str(el.('PolynomB')(ord(1))*factorial(ord(1)-1),format) '; \n'...
 %                   'KNL= ' num2str(factorial(ord(1)-1)*el.('PolynomB')(ord(1)),format) ', \n'...
                %' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL+1),'%d') ' \n'...
                ];
            elelat=[elelat mp '\n'];
        case 'Monitor' % bpm
            pu=[el.('FamName') ' :MONITOR' '; '...\n
                ];
            elelat=[elelat pu '\n'];
        case {'Marker','SkewQuadrupole'}% marker
            mrk=[el.('FamName') ' :MARKER' '; '...\n
                ];
            elelat=[elelat mrk '\n'];
        
       case 'Corrector' % kicker
            ki=[el.('FamName') ' :KICKER' '; '...\n
                ];
            
            elelat=[elelat ki '\n'];
        case 'Drift' % drift
            dr=[el.('FamName') ' : DRIFT, L= ' num2str(el.('Length'),format) '; '];%\n
            elelat=[elelat dr '\n']; 	
         case 'RFCavity' % drift
            rfc=['RF_ON:=1; \n'...
                ...
                el.('FamName') ' : RFCavity, L=' num2str(el.('Length'),format)...
                ',VOLT=' num2str(el.Voltage/1e6,format) ' &\n'... %MV
                ', harm=' num2str(el.HarmNumber,'%3.3d') ', lag=0.0' '; \n'];
      
            elelat=[elelat rfc '\n'];
           
        otherwise
            warning(['Element: ' el.('FamName') ' was not converted, since it does not match any Class.'])
    end
   end
end

def=elelat;

elelat=[];

elelat=[elelat '! LINE \n\n'];

elelat=[elelat linename ' : LINE = ('];

%% define lattice line
% loop all elements
for i=1:length(AT_ring)
    if i==1
        elelat=[elelat '' AT_ring{i}.('FamName') ' &\n'];
    else
        if mod(i,5)==0
        elelat=[elelat ',' AT_ring{i}.('FamName') ' &\n']; 
        else
        elelat=[elelat ',' AT_ring{i}.('FamName') ' ']; 
        end
    end
end

elelat=[elelat ') \n'];

%elelat=[elelat ') \n use ' linename '\n \n twiss, chrom, save \n'];

lines=elelat;

elelat=[def lines];

%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return