function delta_max_rf = atRFacc(ring)
%ATRFACC Computes RF acceptance of the ring
% delta_max_rf = atRFacc(ring)
%   The functions computes the RF acceptance of the ring
%   ring is tha at lattice without radiation
%   delta_max_rf is the RF acceptance
%   
%  See also RFacc

U0=atgetU0(ring);
E0=atenergy(ring);
alpha=mcf(ring,0);
maskcav=atgetcells(ring,'Class','RFCavity');
Voltage=atgetfieldvalues(ring,maskcav,'Voltage');
hnum=atgetfieldvalues(ring,maskcav,'HarmNumber');
h=hnum(1);
Vrf=sum(Voltage);
delta_max_rf=RFacc(Vrf,U0,E0,h,alpha);


