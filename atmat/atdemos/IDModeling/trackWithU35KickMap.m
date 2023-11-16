%read in ESRF lattice
r=atradoff(esrf);

%create element from kick map and plot kick map
U35elem=atidtable_dat('U35',1,'U35.mat',6.04,'IdTablePass');
xtable=U35elem.xtable;
ytable=U35elem.ytable;
xkick2=U35elem.xkick2;
ykick2=U35elem.ykick2;
surf(xtable,ytable,xkick2)
figure
surf(xtable,ytable,ykick2)

%add kick map to lattice
esrf_U35 = add_ID(r,U35elem);

%compute optics with and without kick map
[p0,t0]=atlinopt(r,0,1);
[pU35,tU35]=atlinopt(esrf_U35,0,1);

%write out tune change due to undulator
t0-tU35
