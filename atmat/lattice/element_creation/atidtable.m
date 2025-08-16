function Elem = atidtable(fname, Nslice, filename, Energy, method)
% ATIDTABLE Creates an ID element
%
% FamName	family name
% Nslice	number of slices (1 means the wiggler is represented by a
%           single kick in the center of the device).
% filename	name of file with wiggler tracking tables.
% Energy    Energy of the machine, needed for scaling
% method    tracking function. Defaults to 'IdTablePass'
%
% The tracking table is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% returns assigned structure with class 'KickMap'.
%
% Note: Use atidtable_dat instead, with sort, transpose,
%       and negative first order sign.

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 23-02-2012:  further modifications by B. Nash: reads in only matlab file
% 13-08-2025:  Replaced by atidtable_dat.
%---------------------------------------------------------------------------
    if nargin < 5, method='IdTablePass'; end
    Elem = atidtable_dat(fname,Nslice,filename,Energy,method,...
                        'sort',1, ...
                        'transpose',1, ...
                        'first_order_sign',-1);
end