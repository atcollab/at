function ring=dba()
%create dba lattice

DR_24 = atdrift('DR_24',2.4);
QD = atquadrupole('QD',0.2,-3.2243,'StrMPoleSymplectic4Pass');
DR_03 = atdrift('DR_03',0.3);
QF = atquadrupole('QF',0.2,4.60156,'StrMPoleSymplectic4Pass');
B=atrbend('B',1,0.7853981634,0,'BndMPoleSymplectic4E2Pass');
QDM = atquadrupole('QDM',0.2,-1.95898,'StrMPoleSymplectic4Pass');
SD = atsextupole('SD',0.1,-35.45865);
DR_04 = atdrift('DR_04',0.4);
DR_031309 = atdrift('DR_031309',0.313086069999997);
SF= atsextupole('SF',0.1,13.50225);
QFM = atquadrupole('QFM',0.199999999999999,3.57926,'StrMPoleSymplectic4Pass');
M = atmarker('M');
 
 
cell={...
  DR_24;QD;QD;DR_03;QF;QF;DR_03;QD;QD;DR_031309;...
           B;DR_031309;QDM;QDM;SD;DR_04;SF;QFM;M;QFM;...
        SF;DR_04;SD;QDM;QDM;DR_031309;B;DR_031309;QD;QD;DR_03;...
        QF;QF;DR_03;QD;QD;DR_24;
    };
ring = repmat(cell,4,1);

end
