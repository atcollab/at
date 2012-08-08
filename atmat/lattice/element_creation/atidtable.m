function Elem = atidtable(fname, Nslice, filename, Energy, method)
% atidtable(FAMNAME,Nslice,filename,Energy,method)
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
% returns assigned structure

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 23-02-2012:  further modifications by B. Nash: reads in only matlab file
%---------------------------------------------------------------------------

if nargin < 5, method='IdTablePass'; end
Elem.FamName=fname;
Elem.Nslice=Nslice;
Elem.PassMethod=method;


factor=1/((Energy/0.299792458)^2);
factor1=-1/((Energy/0.299792458));
    
% Read the file, first check if its a matlab file
%[pathstr, name, ext] = fileparts(filename)
%if  ~isequal(ext,'.mat');
 
    D=load(filename);
    xtable=(D.xtable)';
    ytable=(D.ytable)';
    xkick1=factor1*D.xkick1;
    ykick1=factor1*D.ykick1;
    xkick=factor*D.xkick;
    ykick=factor*D.ykick;
    L=D.Len;
% Sort arrays in ascending order and transpose (needed for "IdTablePass.c")

[xtable indx]=sort(xtable); %#ok<TRSRT>
[ytable indy]=sort(ytable); %#ok<TRSRT>

xkick=xkick(indx,indy);
ykick=ykick(indx,indy);

xkick=xkick';
ykick=ykick';
xkick1=xkick1';
ykick1=ykick1';


Elem.Length= L;
Elem.xtable = xtable;
Elem.ytable = ytable;
Elem.xkick = xkick;
Elem.ykick = ykick;
Elem.xkick1 = xkick1;
Elem.ykick1 = ykick1;
