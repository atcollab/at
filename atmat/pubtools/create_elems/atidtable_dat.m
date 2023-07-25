function Elem = atidtable_dat(FamName, Nslice, filename, Energy, method)
% atidtable(FamName, Nslice, filename, Energy, method)
%
% FamName   family name
% Nslice    number of slices (1 means the Insertion Device is represented by a
%           single kick in the center of the device).
% filename  name of file with the ID tracking tables.
% Energy    Normalization energy in GeV, needed for scaling
% method    name of the function to use for tracking. Use 'IdTablePass'
%             or 'WigTablePass'
%
% The tracking table method is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% returns atinsertiondevicekickmap

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 18-07-2023:  blanco-garcia. Eliminates R1,R1,T1,T2,NumIntSteps,MaxOrder
%                               PolynomA, PolynomB,NumX,NumY
%                               (they are not used in IdTablePass)
%                             Adds InsertionDeviceKickMap class
%---------------------------------------------------------------------------

%Elem.FamName        = fname;  % add check for identical family names
%Elem.Filename_in    = filename;
%Elem.Normalization_energy = Energy;
%Elem.Class          = "KickMap";
%Elem.Nslice         = Nslice;
%Elem.MaxOrder			= 3;
%Elem.NumIntSteps 	= 10;
%Elem.R1             = diag(ones(6,1));
%Elem.R2             = diag(ones(6,1));
%Elem.T1             = zeros(1,6);
%Elem.T2             = zeros(1,6);
%Elem.PassMethod 	= method;


factor=1/((Energy/0.299792458)^2);
factor1=-1/((Energy/0.299792458));

% Read the file
D=importdata(filename);
if isfield(D,'Kick1x')
    x=(D.x)';
    y=(D.y)';
    xkick1=factor1*D.Kick1x;
    ykick1=factor1*D.Kick1y;
    xkick=factor*D.Kick2x;
    ykick=factor*D.Kick2y;
    L=D.Len;
    nn=size(xkick1);
    Ny=nn(1);
    Nx=nn(2);
%     ElemData.MultiKick= 1;
%     ElemData.nkicks= nn(3);
else
    [header_mat, data_mat]=mhdrload_bis(filename);
    L=data_mat(:,:,1);
    Nx=data_mat(:,:,2);
    Ny=data_mat(:,:,3);
    A = importdata(filename,' ',10);
    x=A.data;
    x=x(1,1:Nx);
    A = importdata(filename,' ',11);
    txkick=A.data;
    y=txkick(1:Ny,1);
    txkick=txkick(:,2:end);
    A=importdata(filename,' ',11+Ny+3);
    tykick=A.data;
    tykick=tykick(:,2:end);
    A=importdata(filename,' ',11+2*Ny+2*3);
    if isstruct(A)
        txkick1=A.data;
        txkick1=txkick1(:,2:end);
    else
        txkick1=0*txkick;
    end
    A=importdata(filename,' ',11+3*Ny+3*3);
    if isstruct(A)
        tykick1=A.data;
        tykick1=tykick1(:,2:end);
    else
        tykick1=0*tykick;
    end

    xkick=factor*txkick;
    ykick=factor*tykick;

    xkick1=factor1*txkick1;
    ykick1=factor1*tykick1;

    % Sort arrays in ascending order (needed for "IdTablePass.c")
    [y indy]=sort(y);
    [x indx]=sort(x);
    x=x';
    xkick=xkick(indy,indx);
    ykick=ykick(indy,indx);

end

%Deprecated
%Elem.Length= L;
%Elem.NumX = Ny;
%Elem.NumY = Nx;
%Elem.xtable = x;
%Elem.ytable = y;
%Elem.xkick = xkick;
%Elem.ykick = ykick;
%Elem.xkick1 = xkick1;
%Elem.ykick1 = ykick1;
%Elem.PolynomA= [0 0 0 0];
%Elem.PolynomB= [0 0 0 0];

Elem = atinsertiondevicekickmap( FamName, ...
                                 method, ...
                                 filename, ...
                                 Energy, ...
                                 Nslice, ...
                                 L, ...
                                 xkick, ...
                                 ykick, ...
                                 xkick1, ...
                                 ykick1, ...
                                 x, ...
                                 y ...
                               );

