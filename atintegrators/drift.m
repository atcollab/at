function z = drift(fname,L,method)
%DRIFT('FAMILYNAME',Length [m],'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length			length[m]
%		PassMethod		name of the function on disk to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

ElemData.FamName = fname;  % add check for identical family names
ElemData.Length = L;
ElemData.PassMethod=method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;


