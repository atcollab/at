function nuamp(ring,fig,SextFam,scale)
%nuamp(fig,SextFam,scale)
xi0=[0.3644 0.33475];
%esrf=atreadbeta('/Users/boaznash/work/current_projects/AT/s13s20thick.str');
%esrf0=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
%esrf0=atfitchrom(esrf0,xi0,{'S13','S20'},'S19');

%vals = 0:0.001:sqrt(0.01);
%amps = vals.^2;
figure(fig)
clf
hold on
vals = 0:0.001:0.022;
atnuampl(ring,vals);
esrf=scalesext(ring,SextFam,scale);
nuampls(ring,vals);
title(SextFam);
xlabel('hor. amp. (m)')