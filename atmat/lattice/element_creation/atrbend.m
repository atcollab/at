function elem=atrbend(fname,varargin)
%ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
%	creates a rectangular bending magnet element with class 'Bend'
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

[rsrc,L,A,K,method]=decodeatargs({0,0,[],'BendLinearPass'},varargin);
[rsrc,L]=getatarg(rsrc,L,'Length');
[rsrc,A]=getatarg(rsrc,A,'BendingAngle');
[rsrc,K]=getatarg(rsrc,K,'K');
[rsrc,method]=getatarg(rsrc,method,'PassMethod');
[rsrc,PolynomB]=getatarg(rsrc,[0 0],'PolynomB');
if ~isempty(K), PolynomB(2)=K; end
[rsrc,EntranceAngle]=getatarg(rsrc,0.5*A,'EntranceAngle');
[rsrc,ExitAngle]=getatarg(rsrc,0.5*A,'ExitAngle');
elem=atbaselem(fname,method,'Class','Bend','Length',L,...
    'BendingAngle',A,'EntranceAngle',EntranceAngle,'ExitAngle',ExitAngle,...
    'K',PolynomB(2),'PolynomB',PolynomB,rsrc{:});
end
