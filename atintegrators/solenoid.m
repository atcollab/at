function z=solenoid(fname,L,KS,method)

% z=solenoid('FAMILYNAME',Length [m],KS,'METHOD')
%	creates a new family in the FAMLIST - a structure with field
%	FamName			family name
%	Length			length[m]
%	KS              solenoid strength KS [rad/m]
%	PassMethod     name of the function to use for tracking
%
%   function returns assigned address in the FAMLIST that uniquely identifies
%   the family
%
%   Additional structures being set up (initialized to default values within this routine):   
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2           		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit

ElemData.FamName = fname;  % add check for existing identical family names
ElemData.Length = L;
ElemData.K         = KS;
ElemData.MaxOrder = 3;
ElemData.NumIntSteps = 10;
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

