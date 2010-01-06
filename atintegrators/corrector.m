function z = corrector(fname,L,kickangle, method)
%CORRECTOR('FAMILYNAME',LENGTH,ANGLE,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length 			is set to 0 for  marker type 
%		KickAngle       [kickx, kicky] in radians (small) - unis of d(x,y)/ds
%       PassMethod		name of the function on disk to use for tracking
%
% returns assigned index in the FAMLIST that is uniquely identifies
% the family

ElemData.FamName = fname;  % add check for identical family names
ElemData.Length = L;
if length(kickangle) == 2
    ElemData.KickAngle = reshape(kickangle,1,2);
else
    ElemData.KickAngle = [0, 0];
    warning(['Family: ',fname,' - KickAngle is not a 2-vector. Set both components to 0']);
end
    
ElemData.PassMethod = method;

global FAMLIST
z = length(FAMLIST)+1; % number of declared families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

