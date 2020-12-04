%SI units
qe=1.602176620e-19;
partmass = 9.10938356e-31;
E=6.04e9;
me=0.5109989461e6; %rest energy of electron in eV
gamma=E/me;
beta_rel=sqrt(1-1/gamma^2);
clight=3e8;

circ = 844.3907;
current = .005;

wakefact = -qe^2/(partmass*gamma*beta_rel^2*clight^2);
intensity = current*circ/(clight*beta_rel*qe);

nslice = 51;


% Now create wake function
xr = 0.1; %wake extends for 10 cm
table_length = 201;
freqx = 10; %BBR freq of 10 GHz
freqy = 10;
freqz = 10;
qx = 1;
qy = 1;
qz = 1;
rx = .5;
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


%Now generate fast ring and add impedance element

ring=esrf;
indcav=findcells(ring,'Class','RFCavity');
cav=ring(indcav(1));
ring(indcav(:))=[];
ring=[cav;ring];

ring=atsetcavity(ring,8e6,1,992);

[fastring,fastringrad]=atfastring(ring);
fastringBBR=[fastringrad;imp_tab_elem];
