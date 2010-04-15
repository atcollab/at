%esrf(:)
xi0=[0.3644 0.33475];
%esrf=atreadbeta('/Users/boaznash/work/current_projects/AT/esrfdata/s13s20thick.str');
esrf=atreadbeta('s13s20thick.str');
esrf0=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
esrf0=atfitchrom(esrf0,xi0,{'S13','S20'},'S19');

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,1,0);
esrf=scalesext(esrf,'S4',.95);
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
nlchromplot(esrf,-.04,.04,30,16,1,1);
title('S4')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,2,0);
esrf=scalesext(esrf,'S6',.95);
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
nlchromplot(esrf,-.04,.04,30,16,2,1);
title('S6')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,3,0);
esrf=scalesext(esrf,'S13',.95);
%esrf=atfitchrom(esrf,[0.3642 0.3346],'S13','S19');
nlchromplot(esrf,-.04,.04,30,16,3,1);
title('S13 (no xi retune)')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,4,0);
esrf=scalesext(esrf,'S19',.95);
%esrf=atfitchrom(esrf,[0.3642 0.3346],'S13','S19');
nlchromplot(esrf,-.04,.04,30,16,4,1);
title('S19 (no xi retune)')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,5,0);
esrf=scalesext(esrf,'S20',.95);
%esrf=atfitchrom(esrf,[0.3642 0.3346],'S13','S19');
nlchromplot(esrf,-.04,.04,30,16,5,1);
title('S20 (no xi retune)')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,6,0);
esrf=scalesext(esrf,'S22',.95);
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
nlchromplot(esrf,-.04,.04,30,16,6,1);
title('S22')

esrf=esrf0;
nlchromplot(esrf,-.04,.04,30,16,7,0);
esrf=scalesext(esrf,'S24',.95);
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
esrf=atfitchrom(esrf,xi0,{'S13','S20'},'S19');
nlchromplot(esrf,-.04,.04,30,16,7,1);
title('S24')