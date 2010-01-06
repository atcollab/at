function z=quadrupole(fname,L,K,method)
%QUADRUPOLE('FAMILYNAME',Length [m],K,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields%		FamName    
%	FamName			family name
%	Length			length[m]
%	K				K-value of the quadrupole
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2        		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit
%	PassMethod     name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

ElemData.FamName = fname;  % add check for existing identical family names
ElemData.Length = L;
ElemData.K         = K;
ElemData.MaxOrder = 3;
ElemData.NumIntSteps = 10;
ElemData.PolynomA= [0 0 0 0];	 
ElemData.PolynomB= [0 K 0 0];
ElemData.R1 = diag(ones(6,1));
ElemData.R2 = diag(ones(6,1));
ElemData.T1 = zeros(1,6);
ElemData.T2 = zeros(1,6);
ElemData.PassMethod=method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

