function Elem=atmultipole(fname,L,PolynomA,PolynomB,method)

% MULTIPOLE('FAMILYNAME',Length [m],PolynomA,PolynomB,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%	FamName			family name
%	Length			length[m]
%	ElemData.PolynomA= skew [dipole quad sext oct];	 
%	ElemData.PolynomB= normal [dipole quad sext oct]; 
%	PassMethod     name of the function to use for tracking
%
%   internally the additional structure fields are set:
%
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2        		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit4
%
% returns assigned address in the FAMLIST that uniquely identifies
% the family


Elem.FamName = fname;  % add check for identical family names
Elem.Length = L;
Elem.MaxOrder = 3;
Elem.NumIntSteps = 10;
Elem.R1 = diag(ones(6,1));
Elem.R2 = diag(ones(6,1));
Elem.T1 = zeros(1,6);
Elem.T2 = zeros(1,6);
Elem.PolynomA= PolynomA;	 
Elem.PolynomB= PolynomB;
Elem.BendingAngle  	= PolynomB(1);
Elem.PassMethod=method;
