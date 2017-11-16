function elem=atvsbend(fname,varargin)
%ATVSBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
%	creates a sector vertical bending magnet element with class 'VertBend'
%		FAMNAME        	family name
%		LENGTH         	length of the arc for an on-energy particle [m]
%		BENDINGANGLE	total bending angle [rad]
%		K				focusing strength, defaults to 0
%		PASSMETHOD      tracking function, defaults to 'BendLinearPass'
%
%ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,A,K,method]=decodeatargs({0,0,[],'VBndMPoleSymplectic4Pass'},varargin);
[L,rsrc]=getoption(rsrc,'Length',L);
[A,rsrc]=getoption(rsrc,'BendingAngle',A);
[K,rsrc]=getoption(rsrc,'K',K);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[PolynomB,rsrc]=getoption(rsrc,'PolynomB',[0 0]);
[cl,rsrc]=getoption(rsrc,'Class','VertBend');
if ~isempty(K), PolynomB(2)=K; end
[EntranceAngle,rsrc]=getoption(rsrc,'EntranceAngle',0);
[ExitAngle,rsrc]=getoption(rsrc,'ExitAngle',0);
elem=atbaselem(fname,method,'Class',cl,'Length',L,...
    'BendingAngle',A,'EntranceAngle',EntranceAngle,'ExitAngle',ExitAngle,...
    'K',PolynomB(2),'PolynomB',PolynomB,rsrc{:});
end
