function Elem=atrbend(fname,L,A,K,method)
%ATRBEND(FAMNAME,LENGTH,BENDINGANGLE,K,PASSMETHOD)
%	creates a rectangular bending magnet element
%		FAMNAME        	family name
%		LENGTH         	length of the arc for an on-energy particle [m]
%		BENDINGANGLE	total bending angle [rad]
%		K				focusing strength, defaults to 0
%		PASSMETHOD      tracking function, defaults to 'BendLinearPass'
%
%See also: ATDRIFT, ATQUADRUPOLE, ATSEXTUPOLE, ATSBEND
%          ATMULTIPOLE, ATTHINMULTIPOLE, ATMARKER, ATCORRECTOR

if nargin < 5, method='BendLinearPass'; end
if nargin < 4, K=0; end

Elem.FamName=fname;
Elem.Length=L;
Elem.PolynomA=[0 0 0];	 
Elem.PolynomB=[0 K 0]; 
Elem.MaxOrder=2;
Elem.NumIntSteps=10;
Elem.BendingAngle=A;
Elem.EntranceAngle=0.5*A;
Elem.ExitAngle=0.5*A;
Elem.K=K;
Elem.PassMethod=method;
Elem.Class='Bend';
