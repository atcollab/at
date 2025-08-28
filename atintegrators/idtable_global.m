function [z, L] = idtable_global(fname, Nslice, filename, Energy, method)
%IDTABLE Creates a RADIA-Map based element, and adds it to global FAMLIST.
% idtable(fname, Nslice, filename, Energy, method)
%
% FamName	family name
% Nslice	number of slices (1 means the wiggler is represented by a
%           single kick in the center of the device).
% filename	name of file with wiggler tracking tables.
% Energy    Energy of the machine, needed for scaling
% method    name of the function to use for tracking. Use 'WigTablePass'
%
% The tracking table is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
%
% Note: Use atidtable_dat instead, with sort, transpose,
%       and negative first order sign.
%       'KickMap' class has been replaced by 'InsertionDeviceKickMap'.

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z. Marti
% 13-08-2025:  Replaced by atidtable_dat.
%---------------------------------------------------------------------------
    ElemData = atidtable_dat(fname,Nslice,filename,Energy,method,...
                        'sort',1, ...
                        'transpose',1, ...
                        'first_order_sign',-1);
    ElemData.MaxOrder		= 3;
    ElemData.NumIntSteps 	= 10;
    ElemData.R1             = diag(ones(6,1));
    ElemData.R2             = diag(ones(6,1));
    ElemData.T1             = zeros(1,6);
    ElemData.T2             = zeros(1,6);
    ElemData.PolynomA       = [0 0 0 0];
    ElemData.PolynomB       = [0 0 0 0];
    L = ElemData.Length;
    
    global FAMLIST
    z = length(FAMLIST)+1;
    FAMLIST{z}.FamName = fname;
    FAMLIST{z}.NumKids = 0;
    FAMLIST{z}.KidsList= [];
    FAMLIST{z}.ElemData= ElemData;
end