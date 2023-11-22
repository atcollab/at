function Elem = atidtable_dat(FamName, Nslice, filename, Energy, method)
% atidtable_dat Read - RADIA kick maps of 1st and 2nd order in energy
%
% Elem = atidtable_dat(FamName, Nslice, filename, Energy, method)
%
% FamName   family name
% Nslice    number of slices (1 means the Insertion Device is represented
%           by a single kick in the center of the device).
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

%--------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 18-07-2023:  blanco-garcia. Eliminates R1,R1,T1,T2,NumIntSteps,MaxOrder
%                               PolynomA, PolynomB,NumX,NumY
%                               (they are not used in IdTablePass)
%                             Adds InsertionDeviceKickMap class
% 13-11-2023:  blanco-garcia. Renaming variables due to Pull Request
%                             https://github.com/atcollab/at/pull/683
%--------------------------------------------------------------------------

% constants
lightspeed = PhysConstant.speed_of_light_in_vacuum.value;
emassGeV = PhysConstant.electron_mass_energy_equivalent_in_MeV.value/1e3;

% energy scaling for 1st order kick-map
factor1 = -1 / (1e9*sqrt(Energy ^2 - emassGeV^2 )/ lightspeed);
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
else
    isIDfileformat = 3; % text file with only TABs or only SPACES
end

% readout variables according to fileformat
switch isIDfileformat
    case {1} % mat file
        x = (D.x)';
        y = (D.y)';
        xkick1 = factor1 * D.Kick1x;
        ykick1 = factor1 * D.Kick1y;
        xkick2 = factor2 * D.Kick2x;
        ykick2 = factor2 * D.Kick2y;
        L  = D.Len;

    case {2} % new variable names
        x = D.xtable;
        y = D.ytable;
        xkick1 = factor1 * D.xkick1;
        ykick1 = factor1 * D.ykick1;
        xkick2 = factor2 * D.xkick2;
        ykick2 = factor2 * D.ykick2;
        L  = D.Length;

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
        xkick2 = factor2 * txkick2;
        ykick2 = factor2 * tykick2;
        xkick1 = factor1 * txkick1;
        ykick1 = factor1 * tykick1;

        % Sort arrays in ascending order (needed for "IdTablePass.c")
        [y, indy] = sort(y);
        [x, indx] = sort(x);
        x = x';
        xkick2 = xkick2(indy, indx);
        ykick2 = ykick2(indy, indx);
    otherwise
        fprintf('Unsupported fileformat.\n')
end

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
end
