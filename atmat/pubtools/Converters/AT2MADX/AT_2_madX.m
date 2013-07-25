function AT_2_madX(AT_ring,linename)
% this functions converts the AT lattice AT_ring in madX format.
% 
% Beta codes are assumed, please modify to suite your needs.
%
% last update begin 2012

outfile='madXelemdef.elem';
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
    switch el.('BetaCode')
        case 'DI' % dipole
            di=[el.('FamName')   ': '...
                ' RBEND,L=' num2str(el.('Length'),format) ', '...
                ' ANGLE = ' num2str(el.('BendingAngle'),format) ', '...
                ' K1 = ' num2str(el.('PolynomB')(2),format) ', '...
                ' E1= ' num2str(0*el.('EntranceAngle'),format) ', '...
                ' E2= ' num2str(0*el.('ExitAngle'),format) '; '...
               ];
           
            elelat=[elelat di '\n']; %#ok<*AGROW>
        case 'QP' % quadrupole
            
            qp=[el.('FamName')   ': '...
                ' QUADRUPOLE,  L= ' num2str(el.('Length'),format)  ','...
                ' K1= ' num2str(el.('PolynomB')(2),format) ';'...
                ];
            
            elelat=[elelat qp '\n'];
        case 'SX' % sextupole
            sx=[el.('FamName')   ': '...
                ' KSEXT,  L= ' num2str(el.('Length'),format)  ', '...
                ' K2= ' num2str(el.('PolynomB')(3),format) ' ;'...
                ];
            elelat=[elelat sx '\n'];
        case 'MP' % bpm
            % el.('PolynomB')
            ord=find(el.('PolynomB'));
            if isempty(ord)
                ord=length(el.('PolynomB'));
            end
                mp=[el.('FamName')   ': '...
                    ' MULT,  L= ' num2str(el.('Length'),format)  ', '...
                    'ORDER= ' num2str(ord(1)-1) ', &\n' ... 
                   'KNL= ' num2str(el.('PolynomB')(ord(1)),format) ' ;'...
 %                   'KNL= ' num2str(factorial(ord(1)-1)*el.('PolynomB')(ord(1)),format) ', \n'...
                %' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL+1),'%d') ' \n'...
                ];
            elelat=[elelat mp '\n'];
        case 'PU' % bpm
            pu=[el.('FamName') ' :MONITOR' ' ;'...
                ];
            elelat=[elelat pu '\n'];
        case 'MRK'% marker
            mrk=[el.('FamName') ' :MARKER' ' ;'...
                ];
            elelat=[elelat mrk '\n'];
        case 'KI' % kicker
            ki=[el.('FamName') ' :KICKER' ' ;'...
                ];
            
            elelat=[elelat ki '\n'];
        case 'DR' % drift
            dr=[el.('FamName') ' : DRIFT, L= ' num2str(el.('Length'),format) ' ;'];
            elelat=[elelat dr '\n']; 	
         case 'CAV' % drift
            rfc=[el.('FamName') ' : RFCavity, L=' num2str(el.('Length'),format)...
                ',VOLT='  num2str(el.('Voltage'),format) ' ;'...
                ', freq=' num2str(el.('Frequency'),format) '' ' ;'];
      
            elelat=[elelat rfc '\n'];
           
              otherwise
    end
   end
end


elelat=[elelat '! LINE \n\n'];

elelat=[elelat linename ' : LINE = (' ' '];

%% define lattice line
% loop all elements
for i=1:length(AT_ring)
    
 elelat=[elelat ',' AT_ring{i}.('FamName') ' '];

end

elelat=[elelat ') ;'];


%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return