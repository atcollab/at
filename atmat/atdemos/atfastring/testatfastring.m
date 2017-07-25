ring=soleil;
indcav=findcells(ring,'Class','RFCavity');
cav=ring(indcav(1));
ring(indcav(:))=[];
ring=[cav;ring];

ring=atsetcavity(ring,8e6,0,992);

[fastring,fastringrad]=atfastring(ring);
fastringLin=fastring;
fastringLin(3)=[];
m66=fastringLin{2}.M66;
MS=symplectify(m66);
fastringLin{2}.M66=MS;

z0=[0,0,0,0,0.03,0]';
zz0=[1e-5;0;1e-5;0;1e-3;0];

z1Lin=ringpass(fastringLin,zz0,1e7);

z1=ringpass(ring,z0,100);
z1fast=ringpass(fastring,z0,100);
z1fastrad=ringpass(fastringrad,z0,10000);

figure
hold on
plot(z1fastrad(5,:),'-k');