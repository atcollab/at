function z = aperture(fname,varargin)
%APERTURE('FAMILYNAME',limits, method)
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%       Limits         (Xmin, Xmax, Ymin, Ymax) [m]
%		PassMethod		name of the function on disk to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

ElemData=ataperture(fname,varargin{:});

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;


