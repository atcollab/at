function R=ParseAtributesMADX_2_AT(R,elname,var,val)
% determines atribute and sets field in sxs{i} structure AT
%
% created 6-sept-2012

%disp(['checking atribute: ' S8STRING((strfind(S8STRING,'=')+1):end-1)]);

%val=S8STRING((strfind(S8STRING,'=')+1):end);

val0=val;
%S8STRING(1:end)
if ~isempty(val)
    
    val=strrep(val,'=','');
    val=strrep(val,',','');
    val=strrep(val,' ','');
    
    var=strrep(var,'=','');
    var=strrep(var,',','');
    var=strrep(var,' ','');
    
    val=strrep(val,'&','0');
%    val=strrep(val,'E','e');
 %  disp([var '= ' val])
 
%  inde=1:6;
%  [a,b]=meshgrid(inde,inde);
%  mstr=arrayfun(@(i,j)['rm' num2str(i) num2str(j)],a(:),b(:),'un',0);
%  [a,b,c]=meshgrid(inde,inde,inde);
%  tmstr=arrayfun(@(i,j,k)['tm' num2str(i) num2str(k) num2str(j)],a(:),b(:),c(:),'un',0);
  
 mstr={'rm11'
    'rm12'
    'rm13'
    'rm14'
    'rm15'
    'rm16'
    'rm21'
    'rm22'
    'rm23'
    'rm24'
    'rm25'
    'rm26'
    'rm31'
    'rm32'
    'rm33'
    'rm34'
    'rm35'
    'rm36'
    'rm41'
    'rm42'
    'rm43'
    'rm44'
    'rm45'
    'rm46'
    'rm51'
    'rm52'
    'rm53'
    'rm54'
    'rm55'
    'rm56'
    'rm61'
    'rm62'
    'rm63'
    'rm64'
    'rm65'
    'rm66'};

tmstr = {'tm111'
    'tm112'
    'tm113'
    'tm114'
    'tm115'
    'tm116'
    'tm211'
    'tm212'
    'tm213'
    'tm214'
    'tm215'
    'tm216'
    'tm311'
    'tm312'
    'tm313'
    'tm314'
    'tm315'
    'tm316'
    'tm411'
    'tm412'
    'tm413'
    'tm414'
    'tm415'
    'tm416'
    'tm511'
    'tm512'
    'tm513'
    'tm514'
    'tm515'
    'tm516'
    'tm611'
    'tm612'
    'tm613'
    'tm614'
    'tm615'
    'tm616'
    'tm121'
    'tm122'
    'tm123'
    'tm124'
    'tm125'
    'tm126'
    'tm221'
    'tm222'
    'tm223'
    'tm224'
    'tm225'
    'tm226'
    'tm321'
    'tm322'
    'tm323'
    'tm324'
    'tm325'
    'tm326'
    'tm421'
    'tm422'
    'tm423'
    'tm424'
    'tm425'
    'tm426'
    'tm521'
    'tm522'
    'tm523'
    'tm524'
    'tm525'
    'tm526'
    'tm621'
    'tm622'
    'tm623'
    'tm624'
    'tm625'
    'tm626'
    'tm131'
    'tm132'
    'tm133'
    'tm134'
    'tm135'
    'tm136'
    'tm231'
    'tm232'
    'tm233'
    'tm234'
    'tm235'
    'tm236'
    'tm331'
    'tm332'
    'tm333'
    'tm334'
    'tm335'
    'tm336'
    'tm431'
    'tm432'
    'tm433'
    'tm434'
    'tm435'
    'tm436'
    'tm531'
    'tm532'
    'tm533'
    'tm534'
    'tm535'
    'tm536'
    'tm631'
    'tm632'
    'tm633'
    'tm634'
    'tm635'
    'tm636'
    'tm141'
    'tm142'
    'tm143'
    'tm144'
    'tm145'
    'tm146'
    'tm241'
    'tm242'
    'tm243'
    'tm244'
    'tm245'
    'tm246'
    'tm341'
    'tm342'
    'tm343'
    'tm344'
    'tm345'
    'tm346'
    'tm441'
    'tm442'
    'tm443'
    'tm444'
    'tm445'
    'tm446'
    'tm541'
    'tm542'
    'tm543'
    'tm544'
    'tm545'
    'tm546'
    'tm641'
    'tm642'
    'tm643'
    'tm644'
    'tm645'
    'tm646'
    'tm151'
    'tm152'
    'tm153'
    'tm154'
    'tm155'
    'tm156'
    'tm251'
    'tm252'
    'tm253'
    'tm254'
    'tm255'
    'tm256'
    'tm351'
    'tm352'
    'tm353'
    'tm354'
    'tm355'
    'tm356'
    'tm451'
    'tm452'
    'tm453'
    'tm454'
    'tm455'
    'tm456'
    'tm551'
    'tm552'
    'tm553'
    'tm554'
    'tm555'
    'tm556'
    'tm651'
    'tm652'
    'tm653'
    'tm654'
    'tm655'
    'tm656'
    'tm161'
    'tm162'
    'tm163'
    'tm164'
    'tm165'
    'tm166'
    'tm261'
    'tm262'
    'tm263'
    'tm264'
    'tm265'
    'tm266'
    'tm361'
    'tm362'
    'tm363'
    'tm364'
    'tm365'
    'tm366'
    'tm461'
    'tm462'
    'tm463'
    'tm464'
    'tm465'
    'tm466'
    'tm561'
    'tm562'
    'tm563'
    'tm564'
    'tm565'
    'tm566'
    'tm661'
    'tm662'
    'tm663'
    'tm664'
    'tm665'
    'tm666'};
    
    switch var
        case {'l','L'}
            R=[R elname '.(''Length'')=' val '; \n'];
        %disp([var '= ' val])
        case {'angle','ANGLE'}
            R=[R elname '.(''BendingAngle'')=' val '; \n'];
        case {'K1','k1'} % quadrupole
            R=[R elname '.(''PolynomB'')(2)=' val '; \n'];
        case {'K1S',''} % skew quadrupole
            R=[R elname '.(''PolynomA'')(2)=' val '; \n'];
        case {'K2','k2'} % sextupole
            R=[R elname '.(''PolynomB'')(3)=' val '; \n'];
        case {'K3','k3'} % octupole
            R=[R elname '.(''PolynomB'')(4)=' val '; \n'];
        case {'KS','ks'} % solenoid
           R=[R elname '.(''K'')(1)=' val '; \n'];
        case {'E1','e1'} % entrance angle
            R=[R elname '.(''EntranceAngle'')=' val '; \n'];
        case {'E2','e2'} %exit angle
            R=[R elname '.(''ExitAngle'')=' val '; \n'];
            % case 'TILT' % tilt
            %R=[R elname '.(''R1'')=' val '; \n'];
            % case 'APERTURE' % max aperture
            %R=[R elname '.('''')=' val '; \n'];
        case {'HGAP','hgap'} % bending magnet gap
            R=[R elname '.(''FullGap'')=2*' val '; \n'];
        case {'fint','FINT'} % bending magnet entrance fringe
            R=[R elname '.(''FringeInt1'')=' val '; \n'];
        case {'fintx','FINTX'} % bending magnet exit fringe
            R=[R elname '.(''FringeInt2'')=' val '; \n'];
        case {'K1L','k1l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(1)=' val '; \n'];
        case {'K2L','k2l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(2)=' val '; \n'];
        case {'K3L','k3l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(3)=' val '; \n'];
        case {'K4L','k4l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(4)=' val '; \n'];
        case {'K5L','k5l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(5)=' val '; \n'];
        case {'K6L','k6l'} % mulitpole component
            R=[R elname '.(''PolynomB'')(6)=' val '; \n'];
        case {'K1SL','k1sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(1)=' val '; \n'];
        case {'K2SL','k2sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(2)=' val '; \n'];
        case {'K3SL','k3sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(3)=' val '; \n'];
        case {'K4SL','k4sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(4)=' val '; \n'];
        case {'K5SL','k5sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(5)=' val '; \n'];
        case {'K6SL','k6sl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')(6)=' val '; \n'];
        case {'KNL','knl'} % mulitpole component skew
            R=[R elname '.(''PolynomB'')=' val0 '; \n'];
        case {'KSL','ksl'} % mulitpole component skew
            R=[R elname '.(''PolynomA'')=' val0 '; \n'];
        case tmstr % Tijk elements
            if strcmp(var(3),'5')
                var(3)='6';
            elseif strcmp(var(3),'6')
                var(3)='5';
            end
            if strcmp(var(4),'5')
                var(4)='6';
            elseif strcmp(var(4),'6')
                var(4)='5';
            end
            if strcmp(var(5),'5')
                var(5)='6';
            elseif strcmp(var(5),'6')
                var(5)='5';
            end
            R=[R elname '.Tijk(' var(3) ',' var(4) ',' var(5) ')=' strrep(val0,',','') '; \n']; %#ok<*ST2NM>
        case mstr % Mij elements
            if strcmp(var(3),'5')
                var(3)='6';
            elseif strcmp(var(3),'6')
                var(3)='5';
            end
            if strcmp(var(4),'5')
                var(4)='6';
            elseif strcmp(var(4),'6')
                var(4)='5';
            end
            
%             % /(1+delta) for MADX to AT coordinates
%             overdelt=1;
%             if strcmp(var(4),'2') || strcmp(var(4),'4')
%             overdelt=
%             end
            
            R=[R elname '.M66(' var(3) ',' var(4) ')=' strrep(val0,',','') '; \n'];
        case {'VOLT','volt'} %
            R=[R elname '.(''Voltage'')=abs(' val ')*1e6; %% convert to MV \n'];
        case {'FREQ','freq'} %
            R=[R elname '.(''Frequency'')=' val '*1e6; %% convert to MHz  \n'];
        case {'tilt','TILT','Tilt'} % mulitpole component skew
            R=[R elname '.(''Tilt'')(1)=' val '; \n'];
        case {'type','Type','TYPE'} % mulitpole component skew
            R=[R elname '.(''Type'')=' val '; \n'];
        case {'HARMON','harmon'} %
            R=[R elname '.(''HarmNumber'')=' val '; \n'];
            %     case '' %
            %          R=[R elname '.(''')=S8STRING;
            %     case '' %
            %          R=[R elname '.(''')=S8STRING;
        otherwise
            % disp([newel ' : unrecognized element attribute']);
    end
end

%R=[R '\n'];
    
