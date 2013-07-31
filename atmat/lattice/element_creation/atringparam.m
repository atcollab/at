function Elem=atringparam(fname,Energy,per)
%ATRINGPARAM(rname,E0,per)
%	creates a RingParameter Element which should go at the beginning of the ring
%
%FNAME		Family name which may be used as name of Ring
%Energy     Energy of electrons
%PER        Periodicity of the ring (=1 if ring is already expanded)
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE

Elem.FamName=fname;
Elem.Energy=Energy;
Elem.Periodicity=per;
Elem.Length=0;
Elem.PassMethod='IdentityPass';
Elem.Class='RingParam';