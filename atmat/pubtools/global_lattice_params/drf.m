function delta_max_rf=drf(ring,Vrf) %#ok<INUSD>

rp=ringpara(ring(:));
U0=rp.U0;
alpha=rp.alphac;
harm=rp.harm;
E0=rp.E0;

delta_max_rf = sqrt(2*U0/pi/alpha/harm/E0)*sqrt( sqrt((Vrf/U0).^2-1) - acos(U0./Vrf));