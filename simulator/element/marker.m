function z = marker(fname,method)
%MARKER('FAMILYNAME','METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length 			is set to 0 for  marker type 
%		PassMethod		name of the function on disk to use for tracking
%							use 'IdentityPass' for markers that have no action
%							such as reference points
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

ElemData.FamName = fname;  % add check for identical family names
ElemData.Length=0;
ElemData.PassMethod=method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

