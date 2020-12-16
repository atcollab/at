function ring=FODO()
QF=atquadrupole('QF',0.5,1.2,'PassMethod','StrMPoleSymplectic4Pass','Energy',3e9);
QD=atquadrupole('QD',0.5,-1.2,'PassMethod','StrMPoleSymplectic4Pass','Energy',3e9);
Bend=atsbend('Bend',1,2*pi/40,'PassMethod','BndMPoleSymplectic4Pass','Energy',3e9);
RFC=atrfcavity('RFCav');
Dr=atdrift('Dr',0.5);
HalfDr=atdrift('Dr',0.25);
% p1Dr=atdrift('Dr',0.1);
p2Dr=atdrift('Dr',0.2);
% p3Dr=atdrift('Dr',0.3);
SF=atsextupole('SF',0.1,1);
SD=atsextupole('SD',0.1,-1);
cell=[{HalfDr};{Bend};{p2Dr};{SF};{p2Dr};{QF};{Dr};{Bend};{p2Dr};{SD};{p2Dr};{QD};{HalfDr}];

ring=repmat(cell,20,1);
ring=[{RFC};ring];
ring=atsetcavity(ring,5e5,0,100);

end
