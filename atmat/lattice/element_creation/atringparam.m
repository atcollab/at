function Elem=atringparam(rname,Energy,per)
%ATRINGPARAM(rname,E0,per)
%	creates a RingParameter Element which should go at the beginning of the ring
%
%RNAME		name of Ring
%Energy     Energy of electrons
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE

Elem.RingName=rname;
Elem.Energy=Energy;
Elem.Periodicity=per;
Elem.Length=0;
Elem.PassMethod='IdentityPass';
Elem.Class='RingParam';
