function Elem = atidtable_dat(FamName, Nslice, filename, Energy, method)
% atidtable_dat Read - RADIA kick maps of 1st and 2nd order in energy
%
% FamName   family name
% Nslice    number of slices (1 means the Insertion Device is represented by a
%           single kick in the center of the device).
% filename  name of file with the ID tracking tables.
% Energy    Normalization energy in GeV, needed for scaling
% method    name of the function to use for tracking. Use 'IdTablePass'
%             or 'WigTablePass'
%
% returns   atinsertiondevicekickmap
%
%
% This function creates an AT element read from an integrated kickmap file.
%
% The tracking table method is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.

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

lightspeed = PhysConstant.speed_of_light_in_vacuum.value * 1e-9;

% energy scaling for 1st order kick-map
factor1 = -1 / ((Energy / lightspeed));
% energy scaling for 2st order kick-map
factor2 = (factor1) ^ 2;

% check the data separator
fid = fopen(filename);
line_ex = fgetl(fid);
while ischar(line_ex)
line_ex = fgetl(fid);
if line_ex(1) == '#'
    continue
end
flagistab = (line_ex==9);
istab = (sum(flagistab));
if istab
    break
end
end
fclose(fid);
% define data separator
if istab
    datasep = '\t';
else
    datasep = ' ';
end


% Read the file
D = importdata(filename);

if isfield(D, 'Kick1x')
    x = (D.x)';
    y = (D.y)';
    xkick1 = factor1 * D.Kick1x;
    ykick1 = factor1 * D.Kick1y;
    xkick2 = factor2 * D.Kick2x;
    ykick2 = factor2 * D.Kick2y;
    L  = D.Len;
    nn = size(xkick1);
    Ny = nn(1);
    Nx = nn(2);
    %     ElemData.MultiKick= 1;
    %     ElemData.nkicks= nn(3);
else
    A = importdata(filename, datasep, 3);
    L = A.data;
    A = importdata(filename, datasep, 5);
    Nx = A.data;
    A = importdata(filename, datasep, 7);
    Ny = A.data;
    A = importdata(filename, datasep, 10);
    x = A.data;
    x = x(1, 1:Nx);
    A = importdata(filename, datasep, 11);
    txkick = A.data;
    y = txkick(1:Ny, 1);
    txkick = txkick(:, 2:end);
    A = importdata(filename, datasep, 11 + Ny + 3);
    tykick = A.data;
    tykick = tykick(:, 2:end);
    A = importdata(filename, datasep, 11 + 2 * Ny + 2 * 3);

    if isstruct(A)
        txkick1 = A.data;
        txkick1 = txkick1(:, 2:end);
    else
        txkick1 = 0 * txkick;
    end

    A = importdata(filename, datasep, 11 + 3 * Ny + 3 * 3);

    if isstruct(A)
        tykick1 = A.data;
        tykick1 = tykick1(:, 2:end);
    else
        tykick1 = 0 * tykick;
    end

    xkick2 = factor2 * txkick;
    ykick2 = factor2 * tykick;

    xkick1 = factor1 * txkick1;
    ykick1 = factor1 * tykick1;

    % Sort arrays in ascending order (needed for "IdTablePass.c")
    [y, indy] = sort(y);
    [x, indx] = sort(x);
    x = x';
    xkick2 = xkick2(indy, indx);
    ykick2 = ykick2(indy, indx);

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

Elem = atinsertiondevicekickmap( ...
    FamName, ...
    method, ...
    filename, ...
    Energy, ...
    Nslice, ...
    L, ...
    xkick2, ...
    ykick2, ...
    xkick1, ...
    ykick1, ...
    x, ...
    y ...
);
