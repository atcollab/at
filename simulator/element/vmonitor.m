function z = vmonitor(fname,method)
% vmonitor('FAMILYNAME','METHOD')
%       creates a new family in the FAMLIST - a structure with fields
%               FamName                 family name
%               Length                  = 0 for  bpm type
%               PassMethod              name of the function on disk to use fortracking
%                                                       use 'IdentityPass' for bpms that have no action
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
% declare bpms in the lattice file  as BPM = vmonitor('BPM','IdentityPass');


ElemData.FamName = fname;  % add check for identical family names
ElemData.Name='VBPM';
ElemData.Length=0;
ElemData.PassMethod=method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;
