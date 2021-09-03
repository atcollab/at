function AT_2_Elegant(AT_ring,linename)
% this functions converts the AT lattice AT_ring in elegant form.
%
% FringeQuadEntranceExit are ignored.
% apertures defined in elements are ignored.
% EDRIFT is used instead of DRIFT. 
% if radiation is on in dipoles, cavity phase is set as in atsetcavity


outfile='elegantconvertedlattice.lte';
elelat=['!!\n!!  Elegant lattice: ' linename '\n!!  Created: ' datestr(now) '\n!!\n!!\n\n'];
nkickperL=40; % kicks per meter

%% get family names for definitions
[families,ind_first_oc_ring]=...
    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');

elelat=[elelat '! DEFINITIONS \n\n'];

form='%8.16f';
C = findspos(AT_ring,length(AT_ring)+1);
[~,~,V_rf,~,U0] = atenergy(AT_ring);


% get passmethods to see if radiation is on.
Radiation = false;
pm = unique(atgetfieldvalues(AT_ring,atgetcells(AT_ring,'PassMethod'),'PassMethod'));
if strfind(pm{cellfun(@(a)not(~contains(a,'Bnd')),pm,'un',1)},'Rad')
   Radiation = true;
end
            
%% loop families for definitions
for i=1:length(families)
    el= AT_ring{ind_first_oc_ring(i)};
    switch el.('Class')
        case 'Bend' % dipole
            
            di=[el.('FamName')   ': &\n'...
                ' CSBEND,L=' num2str(el.('Length'),form) ', &\n'...
                ...   %             ' CSBEND,L=' num2str((el.('Length')/2)*el.('BendingAngle')/(sin(el.('BendingAngle')/2)),form) ', &\n'...
                ' ANGLE = ' num2str(el.('BendingAngle'),form) ', &\n'...
                ' K1 = ' num2str(el.('PolynomB')(2),form) ', &\n'...
                ' SYNCH_RAD= ' num2str(contains(el.('PassMethod'),'Rad'),'%d') ', &\n'...
                ' E1= ' num2str(el.('EntranceAngle'),form) ', &\n'...
                ' E2= ' num2str(el.('ExitAngle'),form) ', &\n' ];
            ... %                ' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL),'%d') ' \n'...
                if(isfield(el,'T1'))
                di= [ di ' DX= ' num2str(-1*el.('T1')(1),form) ', &\n'...
                    ' DY= ' num2str(-1*el.('T1')(3),form) ', &\n' ] ;
                end
                di= [di ...
                    ' N_KICKS= 40 \n'...
                    ];
                elelat=[elelat di '\n']; %#ok<*AGROW>
                
        case 'Corrector' % Hor Corr
            hc= [el.('FamName') ' : KICKER \n'];
            elelat=[elelat hc '\n'];
            
        case 'Quadrupole' % quadrupole
            
            qp=[el.('FamName')   ': &\n'...
                ' KQUAD,  L= ' num2str(el.('Length'),form)  ', &\n'...
                ' SYNCH_RAD= ' num2str(contains(el.('PassMethod'),'Rad'),'%d') ', &\n'...
                ' K1= ' num2str(el.('PolynomB')(2),form) ', &\n'];
            ... %                 ' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL),'%d') ' \n'...
                if(isfield(el,'Tilt'))
                qp= [qp ' TILT= ' num2str(el.('Tilt')) ' & \n'] ;
                end
                
                if(isfield(el,'T1'))
                    qp = [ qp ' DX= ' num2str(-1*el.('T1')(1),form) ', &\n'...
                        ' DY= ' num2str(-1*el.('T1')(3),form) ', &\n' ] ;
                end
                qp= [qp ...
                    ' N_KICKS= 40 \n'...
                    ];
                
                elelat=[elelat qp '\n'];
        case 'Sextupole' % sextupole
            sx=[el.('FamName')   ': &\n'...
                ' KSEXT,  L= ' num2str(el.('Length'),form)  ', &\n'...
                ...%                ' K2= ' num2str(el.('PolynomB')(3),form) ', &\n'...
                ' SYNCH_RAD= ' num2str(contains(el.('PassMethod'),'Rad'),'%d') ', &\n'...
                ' K2= ' num2str(2*el.('PolynomB')(3),form) ', &\n'];
            ... %                ' N_KICKS= ' num2str(ceil(el.('Length')*nkickperL),'%d') ' \n'...
                if(isfield(el,'T1'))
                sx=[sx ' DX= ' num2str(-1*el.('T1')(1),form) ', &\n'...
                    ' DY= ' num2str(-1*el.('T1')(3),form) ', &\n'];
                end
                if(isfield(el,'Tilt'))
                    sx = [ sx ' TILT= ' num2str(el.('Tilt')) ', &\n' ];
                end
                sx = [ sx ' N_KICKS= 40 \n'];
                
                elelat=[elelat sx '\n'];
        case 'Octupole' % octupole
            sx=[el.('FamName')   ': &\n'...
                ' KOCT,  L= ' num2str(el.('Length'),form)  ', &\n'...
                ' SYNCH_RAD= ' num2str(contains(el.('PassMethod'),'Rad'),'%d') ', &\n'...
                ...%                ' K2= ' num2str(el.('PolynomB')(3),form) ', &\n'...
                ' K3= ' num2str(6*el.('PolynomB')(4),form) ', &\n'];
            if(isfield(el,'T1'))
                sx=[sx ' DX= ' num2str(-1*el.('T1')(1),form) ', &\n'...
                    ' DY= ' num2str(-1*el.('T1')(3),form) ', &\n'];
            end
            if(isfield(el,'Tilt'))
                sx = [ sx ' TILT= ' num2str(el.('Tilt')) ', &\n' ];
            end
            sx = [ sx ' N_KICKS= 40 \n'];
            
            elelat=[elelat sx '\n'];
        case 'Multipole'
            % el.('PolynomB')
            ord=find(el.('PolynomB'));
            if isempty(ord)
                ord=length(el.('PolynomB'));
            end
            mp=[el.('FamName')   ': &\n'...
                ' MULT,  L= ' num2str(el.('Length'),form)  ', &\n'...
                'ORDER= ' num2str(ord(1)-1) ', &\n' ...
                'SYNCH_RAD= ' num2str(contains(el.('PassMethod'),'Rad'),'%d') ', &\n'...
                'KNL= ' num2str(factorial(ord-1)*el.('Length')*el.('PolynomB')(ord(1)),form) ', &\n'];
            ...%                   'KNL= ' num2str(factorial(ord(1)-1)*el.('PolynomB')(ord(1)),form) ', \n'...
                if(isfield(el,'T1'))
                mp=[mp ' DX= ' num2str(-1*el.('T1')(1),form) ', &\n'...
                    ' DY= ' num2str(-1*el.('T1')(3),form) ', &\n'];
                end
                if(isfield(el,'Tilt'))
                    mp = [ mp ' TILT= ' num2str(el.('Tilt')) ', &\n' ];
                end
                mp = [ mp ' N_KICKS= 40 \n'];
                
                elelat=[elelat mp '\n'];
        case 'Monitor' % bpm
            pu=[el.('FamName') ' : MONI, L= ' num2str(el.('Length'),form) ' \n'];
            elelat=[elelat pu '\n'];
        case 'Marker'% marker
            mrk=[el.('FamName') ' :MARK' ' \n'...
                ];
            elelat=[elelat mrk '\n'];
        case 'Kicker' % kicker
            ki=[el.('FamName') ' :KICKER' ' \n'...
                ];
            
            elelat=[elelat ki '\n'];
        case 'Drift' % drift
            dr=[el.('FamName') ' : EDRIFT, L= ' num2str(el.('Length'),form) ' \n'];
            elelat=[elelat dr '\n'];
        case 'RFCavity' % drift
            if Radiation
            rfc=[el.('FamName') ' : RFCA, L=' num2str(el.('Length'),'%8.12f')...
                ',VOLT='  num2str(el.('Voltage'),'%8.6f') ' &\n'...
                ',phase="180 ' num2str(U0) ' ' num2str(V_rf) ' / dasin -",' ' &\n'...
                ' freq="c_mks ' num2str(C,form) ' / ' num2str(el.('HarmNumber'),'%d') ' *" ' ' \n'];
            else
            rfc=[el.('FamName') ' : RFCA, L=' num2str(el.('Length'),'%8.12f')...
                ',VOLT='  num2str(el.('Voltage'),'%8.6f') ' &\n'...
                ',phase="180",' ' &\n'...
                ' freq="c_mks ' num2str(C,form) ' / ' num2str(el.('HarmNumber'),'%d') ' *" ' ' \n'];
            end
            
            elelat=[elelat rfc '\n'];
            
        otherwise
    end
    
end

otherlines=['!Malign and watch \n'...
    'M1: MALIGN,on_pass=0\n'...
    'W1: watch, filename="%%s.w1", mode ="centroid"\n\n'...
    '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!  \n'...
    '!	full 32-cell ring(line):			! \n'...
    '!	             						! \n'...
    '!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! \n\n'
    ];
elelat=[elelat otherlines '\n'];


elelat=[elelat '! LINE \n\n'];

elelat=[elelat linename ' : LINE = (M1,W1,' ' &\n'];

%% define lattice line
% loop all elements
for i=1:(length(AT_ring)-1)
    
    elelat=[elelat AT_ring{i}.('FamName') ', '];
    if (floor(i/8)==ceil(i/8))
        elelat = [elelat '&\n'];
    end
end

elelat=[elelat AT_ring{length(AT_ring)}.('FamName')  ') \n\n'];


%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return