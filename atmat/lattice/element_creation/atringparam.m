function elem=atringparam(fname,varargin)
%ATRINGPARAM(FAMNAME,E0,NBPERIODS)
%	creates a RingParameter Element which should go at the beginning of the ring
%
%FAMNAME	Family name which may be used as name of Ring
%E0         Energy of electrons
%NBPERIODS	Periodicity of the ring (1 if ring is already expanded)
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE

[rsrc,energy,nbper]=decodeatargs({6E9,1},varargin);
[rsrc,energy]=getatarg(rsrc,energy,'Energy');
[rsrc,nbper]=getatarg(rsrc,nbper,'Periodicity');
elem=atbaselem(fname,'IdentityPass','Class','RingParam',...
    'Length',0,...
    'Energy',energy,'Periodicity',nbper,rsrc{:});
end
