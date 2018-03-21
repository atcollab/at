function elem=atringparam(fname,varargin)
%ATRINGPARAM Creates a RingParameter Element which should go at the beginning of the ring
%
%  atringparam(FAMNAME,E0,NBPERIODS)
%	
%  INPUTS
%  1. FAMNAME	- Family name which may be used as name of Ring
%  2. E0        - Energy of electrons
%  3. NBPERIODS - Periodicity of the ring (1 if ring is already expanded)
%
%  OUTPUTS
%  1. elem - RingParam class elem
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atrbend
%          atmultipole, atthinmultipole

[rsrc,energy,nbper]=decodeatargs({6E9,1},varargin);
[energy,rsrc]=getoption(rsrc,'Energy',energy);
[nbper,rsrc]=getoption(rsrc,'Periodicity',nbper);
[cl,rsrc]=getoption(rsrc,'Class','RingParam');
elem=atbaselem(fname,'IdentityPass','Class',cl,...
    'Energy',energy,'Periodicity',nbper,rsrc{:});
end
