%read in ESRF lattice
r=atradoff(esrf);

% create element from kick map and plot kick map.
% sort and transpose are needed for U35.mat
U35elem=atidtable_dat('U35',1,'U35.mat',6.04,'IdTablePass', ...
                      'sort',1,'transpose',1);
xtable=U35elem.xtable;
ytable=U35elem.ytable;
xkick2=U35elem.xkick;
ykick2=U35elem.ykick;
surf(xtable,ytable,xkick2)
figure
surf(xtable,ytable,ykick2)

%add kick map to lattice
esrf_U35 = add_ID(r,U35elem);

%compute optics with and without kick map
[p0,t0]=atlinopt(r,0,1);
[pU35,tU35]=atlinopt(esrf_U35,0,1);

%write out tune change due to undulator
deltatune = t0-tU35;
fprintf("delta_tune = %.3g, %.3g\n",deltatune(1),deltatune(2));
