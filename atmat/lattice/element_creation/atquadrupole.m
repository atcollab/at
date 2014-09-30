function elem=atquadrupole(fname,varargin)
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD)
%	creates a quadrupole element with Class 'Quadrupole'
%
%FAMNAME	family name
%LENGTH     length [m]
%K          strength [m-2]
%PASSMETHOD	tracking function, defaults to 'QuadLinearPass'
%
%ATQUADRUPOLE(FAMNAME,LENGTH,K,PASSMETHOD,'FIELDNAME1',VALUE1,...)
%   Each pair {'FIELDNAME',VALUE} is added to the element
%
%See also: ATDRIFT, ATSEXTUPOLE, ATSBEND, ATRBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

[rsrc,L,K,method]=decodeatargs({0,[],'QuadLinearPass'},varargin);
[rsrc,L]=getatarg(rsrc,L,'Length');
[rsrc,K]=getatarg(rsrc,K,'K');
[rsrc,method]=getatarg(rsrc,method,'PassMethod');
[rsrc,PolynomB]=getatarg(rsrc,[0 0],'PolynomB');
[rsrc,maxorder]=getatarg(rsrc,1,'MaxOrder');
if ~isempty(K), PolynomB(2)=K; end
elem=atbaselem(fname,method,'Class','Quadrupole','Length',L,...
    'K',PolynomB(2),'PolynomB',PolynomB,'MaxOrder',maxorder,rsrc{:});
end
