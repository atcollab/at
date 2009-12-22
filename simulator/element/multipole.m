function z=multipole(fname,L,PolynomA,PolynomB,method)

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


ElemData.FamName = fname;  % add check for identical family names
ElemData.Length = L;
ElemData.MaxOrder = 3;
ElemData.NumIntSteps = 10;
ElemData.R1 = diag(ones(6,1));
ElemData.R2 = diag(ones(6,1));
ElemData.T1 = zeros(1,6);
ElemData.T2 = zeros(1,6);
ElemData.PolynomA= PolynomA;	 
ElemData.PolynomB= PolynomB;
ElemData.BendingAngle  	= PolynomB(1);
ElemData.PassMethod=method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

