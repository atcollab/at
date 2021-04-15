function imp_tab_elem=makeBBR(cicumference, bunch_current)
%SI units;

clight = PhysConstant.speed_of_light_in_vacuum.value;   % m/s
qe=PhysConstant.elementary_charge.value;                % C
partmass = PhysConstant.electron_mass.value;            % kg
me = 1.0e6 * PhysConstant.electron_mass_energy_equivalent_in_MeV.value; % eV
E = 6.04e9;                                             % eV
gamma=E/me;
beta_rel=sqrt(1-1/gamma^2);

wakefact = -qe^2/(partmass*gamma*beta_rel^2*clight^2);
intensity = bunch_current*cicumference/(clight*beta_rel*qe);

nslice = 51;

% Now create wake function
xr = 0.1;       % wake extends for 10 cm
table_length = 201;
freqx = 10;     % BBR freq of 10 GHz
freqy = 10;
freqz = 10;
qx = 1;
qy = 1;
qz = 1;
rx = .5;        % MOhm
ry = 2;
rz = .01;

[ s,bbrx,bbry,bbrz ] = bbr_gentab(xr,table_length,freqx,freqy,freqz,qx,qy,qz,rx,ry,rz);

betax_obs = 1; %get from lattice
betay_obs = 1;

imp_tab_elem=atbaselem('imp_tab','ImpedanceTablePass');
imp_tab_elem.Nslice = nslice; %101
imp_tab_elem.Intensity = intensity;% (number of charges/bunch 7.0e10 - MW_th @ esrf)
imp_tab_elem.Wakefact = wakefact;
imp_tab_elem.Nelem = table_length;
imp_tab_elem.WakeT = s;
imp_tab_elem.WakeDX = bbrx; %[V/C/m]
imp_tab_elem.WakeDY = bbry; %[V/C/m]
imp_tab_elem.WakeQX = zeros(table_length,1);
imp_tab_elem.WakeQY = zeros(table_length,1);
imp_tab_elem.WakeZ = bbrz; %[V/C]
imp_tab_elem.On_x = 0.0;
imp_tab_elem.On_y = 0.0;
imp_tab_elem.On_z = 1.0;
imp_tab_elem.On_qx = 0.0;
imp_tab_elem.On_qy = 0.0;
imp_tab_elem.Normfactx=1.0/betax_obs;
imp_tab_elem.Normfacty=1.0/betay_obs;

    function [ s,bbrx,bbry,bbrz ] = bbr_gentab(xr,npoints,freqx,freqy,freqz,qx,qy,qz,rx,ry,rz)
        %xr is extent of wake function in m
        %npoints is number of points in table
        %freq in GHz
        %Rshunt in MOhm
        
        %Based on equations 2.84 and 2.88 in Chao's 'Physics of Collective Instabilities'
        
        clt = PhysConstant.speed_of_light_in_vacuum.value;
        if mod(npoints,2) == 0
            npoints=npoints+1;
        end
        
        s = linspace(-xr,xr,npoints);
        
        omegarx = 2.0*pi*freqx*1.0e9;
        alphax = omegarx/(2.0*qx);
        omegabarx = sqrt(omegarx*omegarx-alphax*alphax);
        
        omegary = 2.0*pi*freqy*1.0e9;
        alphay = omegary/(2.0*qy);
        omegabary = sqrt(omegary*omegary-alphay*alphay);
        
        omegarz = 2.0*pi*freqz*1.0e9;
        alphaz = omegarz/(2.0*qz);
        omegabarz = sqrt(omegarz*omegarz-alphaz*alphaz);
        
        bbrx=zeros(length(s),1);
        bbry=zeros(length(s),1);
        bbrz=zeros(length(s),1);
        
        for i = 1:length(s)
            %chao 2.88
            %chao 2.84
            if s(i)==0.0
                bbrz(i)=rz*1e6*alphaz;
                bbrx(i)=0.0;
                bbry(i)=0.0;
            elseif s(i)<0.0
                bbrx(i)=0.0;
                bbry(i)=0.0;
                bbrz(i)=0.0;
            else
                bbrx(i) = rx*1e6*omegarx*omegarx/qx/omegabarx*exp(alphax*(-s(i))/clt)*sin(omegabarx*(-s(i))/clt);
                bbry(i) = ry*1e6*omegary*omegary/qy/omegabary*exp(alphay*(-s(i))/clt)*sin(omegabary*(-s(i))/clt);
                bbrz(i) = rz*1e6*2.0*alphaz*exp(alphaz*(-s(i))/clt)*(cos(omegabarz*(-s(i))/clt) + ...
                    alphaz/omegabarz*sin(omegabarz*(-s(i))/clt));
            end
        end
        
    end
end
