function Elem = atidtable_dat(FamName, Nslice, filename, Energy, method,...
                              varargin)
% atidtable_dat Creates an element from RADIA files containing kick maps
%               of second order in energy, and optionally first order.
%
% Elem = atidtable_dat(FamName, Nslice, filename, Energy, method)
%
% Arguments:
% FamName   family name
% Nslice    number of slices (1 means the Insertion Device is represented
%           by a single kick in the center of the device).
% filename  name of file with the ID tracking tables.
% Energy    Normalization energy in GeV, needed for scaling
% method    name of the function to use for tracking. Use 'IdTablePass'.
%
% Options:
% sort      sort the imput table if not zero. Default, 0.
% transpose transpose the imput table if not zero. Default, 0.
% verbose   print info if not zero. Default, 0.
% first_order_sign Defines the sign of first order table. Default, +1.
%
% Returns
%   atinsertiondevicekickmap
%
% This elememt implements tracking through an integrated magnetic
% field map of second order in energy, normalized to and energy
% value that is required to calculate alpha. See Eq. (5) in [#].
%
% First order maps could be included. See Eq. (3) in [#].
% Note that positive and negative signs are not taken into account in
% this implementation. Input should already include the sign difference.
%
% Default PassMethod: ``IdTablePass``.
%
% [#] Pascale ELLEAUME, "A New Approach to the Electron Beam Dynamics in
%     Undulators and  Wigglers". EPAC1992 0661.
%     European Synchrotron Radiation Facility.
%     BP 220, F-38043 Grenoble, France

%--------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 18-07-2023:  oblanco. Eliminates R1,R1,T1,T2,NumIntSteps,MaxOrder
%                               PolynomA, PolynomB,NumX,NumY
%                               (they are not used in IdTablePass)
%                             Adds InsertionDeviceKickMap class
% 15-08-2025   oblanco. Eliminates redundant idtable_global.m and
%                           atidtable.m.
%                       Add options for compatibility.
%--------------------------------------------------------------------------

p = inputParser;
addOptional(p,'sort',0);
addOptional(p,'transpose',0);
addOptional(p,'verbose',0);
addOptional(p,'first_order_sign',1);
parse(p,varargin{:});
par = p.Results;

dosort = par.sort;
dotranspose = par.transpose;
verbose = par.verbose;
first_order_sign = par.first_order_sign;

% constants
lightspeed = PhysConstant.speed_of_light_in_vacuum.value;
emassGeV = PhysConstant.electron_mass_energy_equivalent_in_MeV.value/1e3;

% energy scaling for 1st order kick-map
factor1 = first_order_sign * ...
            1 / (1e9*sqrt(Energy ^2 - emassGeV^2 )/ lightspeed);
% energy scaling for 2st order kick-map
factor2 = (factor1) ^ 2;

% Read the file
D = importdata(filename);

% mat file ?
if isfield(D, 'Kick1x') % old matlab format
    isIDfileformat = 1;
    fprintf('Old format type. Variables will be renamed.\n');
elseif isfield(D, 'xkick1') % ID new variable names
    isIDfileformat = 2;
    if verbose, fprintf('Found xkick1 in .mat file.\n'); end
else
    isIDfileformat = 3; % text file with only TABs or only SPACES
    if verbose, fprintf('Text file.\n'); end
end

% readout variables according to fileformat
switch isIDfileformat
    case {1} % mat file
        x = (D.x)';
        y = (D.y)';
        xkick1 = factor1 * D.Kick1x;
        ykick1 = factor1 * D.Kick1y;
        xkick = factor2 * D.Kick2x;
        ykick = factor2 * D.Kick2y;
        L  = D.Len;

    case {2} % new variable names
        x = D.xtable;
        y = D.ytable;
        xkick1 = factor1 * D.xkick1;
        ykick1 = factor1 * D.ykick1;
        xkick = factor2 * D.xkick;
        ykick = factor2 * D.ykick;
        if isfield(D,'Length'), L  = D.Length; end
        if isfield(D,'Len'), L  = D.Len; end

    case {3}% data from text file
        % check the data separator, either all tabs or all spaces
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

        % read data
        A = importdata(filename, datasep, 3);
        L = A.data;
        A = importdata(filename, datasep, 5);
        Nx = A.data;
        A = importdata(filename, datasep, 7);
        Ny = A.data;
        % kick tables
        A = importdata(filename, datasep, 10);
        x = A.data;
        x = x(1, 1:Nx);
        A = importdata(filename, datasep, 11);
        txkick2 = A.data;
        y = txkick2(1:Ny, 1);
        txkick2 = txkick2(:, 2:end);
        A = importdata(filename, datasep, 11 + Ny + 3);
        tykick2 = A.data;
        tykick2 = tykick2(:, 2:end);

        % check if first orders kicks are defined, otherwise set them to zero
        A = importdata(filename, datasep, 11 + 2 * Ny + 2 * 3);
        if isstruct(A)
            txkick1 = A.data;
            txkick1 = txkick1(:, 2:end);
        else
            txkick1 = 0 * txkick2;
        end
        A = importdata(filename, datasep, 11 + 3 * Ny + 3 * 3);
        if isstruct(A)
            tykick1 = A.data;
            tykick1 = tykick1(:, 2:end);
        else
            tykick1 = 0 * tykick2;
        end

        % scale kick tables
        xkick = factor2 * txkick2;
        ykick = factor2 * tykick2;
        xkick1 = factor1 * txkick1;
        ykick1 = factor1 * tykick1;

        % Sort arrays in ascending order (needed for "IdTablePass.c")
        [y, indy] = sort(y);
        [x, indx] = sort(x);
        x = x';
        xkick = xkick(indy, indx);
        ykick = ykick(indy, indx);
    otherwise
        fprintf('Unsupported fileformat.\n');
        return;
end

if dosort
    if verbose, fprintf('Sorting data.\n'); end
    [x, indx]=sort(x);
    [y, indy]=sort(y);
    xkick=xkick(indx,indy);
    ykick=ykick(indx,indy);
end

if dotranspose
    if verbose, fprintf('Transposing data.\n'); end
    xkick=xkick';
    ykick=ykick';
    xkick1=xkick1';
    ykick1=ykick1';
end

Elem = atinsertiondevicekickmap( ...
    FamName, ...
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
end
