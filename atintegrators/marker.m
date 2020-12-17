function z = marker(fname,method)
%MARKER Creates a marker element in old AT version (Obsolete)
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
%
%  NOTES
%  1. Obsolete: use atmarker instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

ElemData = atmarker(fname,method);

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

