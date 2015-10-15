esrf=build_simple(1,'s13s20thick.str');
indcav=findcells(esrf,'Class','RFCavity');
cav=esrf(indcav(1));
esrf(indcav(:))=[];
esrf=[cav;esrf];

esrf=atsetcavity(esrf,8e6,0,992);

[fastring,fastringrad]=atfastring2(esrf);
fastringLin=fastring;
fastringLin(3)=[];
m66=fastringLin{2}.M66;
MS=symplectify(m66);
fastringLin{2}.M66=MS;

z0=[0,0,0,0,0.03,0]';
zz0=[1e-5;0;1e-5;0;1e-3;0];

z1Lin=ringpass(fastringLin,zz0,1e7);

z1=ringpass(esrf,z0,100);
z1fast=ringpass(fastring,z0,100);
z1fastrad=ringpass(fastringrad,z0,10000);

figure
hold on
%plot(z1(1,:),z1(2,:),'.r');
%plot(z1fast(1,:),z1fast(2,:),'.b');
%plot(z1fastrad(1,:),z1fastrad(2,:),'.k');
plot(z1fastrad(5,:),'-k');