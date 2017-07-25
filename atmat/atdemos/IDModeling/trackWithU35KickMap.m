%read in ESRF lattice
r=esrf;

%create element from kick map and plot kick map
U35elem=atidtable('U35',1,'U35.mat',6.04,'IDTablePass');
xtable=U35elem.xtable;
ytable=U35elem.ytable;
xkick=U35elem.xkick;
ykick=U35elem.ykick;
surf(xtable,ytable,xkick)
figure
surf(xtable,ytable,ykick)

%add kick map to lattice
esrf_U35 = add_ID(r,U35elem);

%compute optics with and without kick map
[p0,t0]=atlinopt(r,0,1);
[pU35,tU35]=atlinopt(esrf_U35,0,1);

%write out tune change due to undulator
t0-tU35