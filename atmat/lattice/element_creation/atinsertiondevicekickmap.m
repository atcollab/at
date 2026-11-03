function Elem = atinsertiondevicekickmap(fname,varargin)
%ATINSERTIONDEVICEKICKMAP Create an insertion device kick-map element
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME)
%   Create an empty zero-length insertion device normalized at 0 GeV, with
%   PASSMETHOD='DriftPass'.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,NSLICE,FILENAME,...
%   NORMALIZATION_ENERGY,[PASSMETHOD],'FIELD1',VALUE1,...)
%   Create an element using the same positional metadata order as pyAT.
%   PASSMETHOD defaults to 'DriftPass'. The kick tables and all other
%   element attributes are supplied as field/value pairs.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,[PASSMETHOD],...
%   'FIELD1',VALUE1,...)
%   Reconstruct an element entirely from field/value pairs, as used by
%   lattice serialization.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,PASSMETHOD,FILENAME,...
%   NORMALIZATION_ENERGY,NSLICE,LENGTH,XKICK,YKICK,XKICK1,YKICK1,...
%   XTABLE,YTABLE)
%   Read the historical full-positional lattice format.
%
% Inputs:
%   FAMNAME              Family name
%   NSLICE               Number of integration slices
%   FILENAME             Source kick-map file name
%   NORMALIZATION_ENERGY Energy in GeV used to normalize the kick map
%   PASSMETHOD           Tracking function. Default: 'DriftPass'
%
% Field/value pairs:
%   'Length'             Insertion device length [m]
%   'xkick', 'ykick'     Second-order horizontal and vertical kick tables
%   'xkick1', 'ykick1'   First-order horizontal and vertical kick tables
%   'xtable', 'ytable'   Horizontal and vertical table coordinates
%   'KickmapStore'       Structure of named kick maps. Each entry contains
%                        Filename_in, Normalization_energy, Nslice, Length,
%                        xkick, ykick, xkick1, ykick1, xtable and ytable
%   'ActiveKickmap'      Name of the KickmapStore entry currently copied to
%                        the element tracking fields
%
% Additional standard AT element fields may also be supplied as field/value
% pairs. PassMethod may equivalently be supplied as a 'PassMethod' field.
%
% Examples:
%   emptyid = atinsertiondevicekickmap('ID');
%   id = atinsertiondevicekickmap('ID',25,'id_kicks.txt',2.75,...
%       'IdTablePass','Length',2.0,...
%       'xkick',xkick,'ykick',ykick,...
%       'xkick1',xkick1,'ykick1',ykick1,'xtable',x,'ytable',y);
%
% The tracking method is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% Returns an element with Class 'InsertionDeviceKickMap'

%------------------
% Modification Log:
% -----------------
% 24-05-2023:  blanco-garcia, added for compatibility with pyat
%---------------------------------------------------------------------------

nslice = 1;
kickmapFilename = '';
normalizationEnergy = 0.0;
args = varargin;

isLegacy = length(args) == 11 && ...
    (ischar(args{1}) || (isstring(args{1}) && isscalar(args{1}))) && ...
    endsWith(args{1},'Pass') && ...
    (ischar(args{2}) || (isstring(args{2}) && isscalar(args{2}))) && ...
    isnumeric(args{3}) && isnumeric(args{4});

if isLegacy
    method = args{1};
    kickmapFilename = args{2};
    normalizationEnergy = args{3};
    nslice = args{4};
    rsrc = {'Length',args{5}, ...
        'xkick',args{6},'ykick',args{7}, ...
        'xkick1',args{8},'ykick1',args{9}, ...
        'xtable',args{10},'ytable',args{11}};
else
    if ~isempty(args) && isnumeric(args{1})
        if length(args) < 3
            error('AT:WrongArgument', ...
                ['NSLICE, FILENAME and NORMALIZATION_ENERGY must be ' ...
                 'supplied together']);
        end
        nslice = args{1};
        kickmapFilename = args{2};
        normalizationEnergy = args{3};
        args(1:3) = [];
    end
    [rsrc,method] = decodeatargs({'DriftPass'},args);
end

[kickmapFilename,rsrc] = getoption(rsrc,'Filename_in',kickmapFilename);
[normalizationEnergy,rsrc] = getoption(rsrc,'Normalization_energy',normalizationEnergy);
[nslice,rsrc] = getoption(rsrc,'Nslice',nslice);
Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap', ...
    'Filename_in',kickmapFilename, ...
    'Normalization_energy',normalizationEnergy, ...
    'Nslice',nslice,rsrc{:});

if ~isfield(Elem,'KickmapStore') && isfield(Elem,'xkick')
    kickmapData = struct( ...
        'Filename_in',Elem.Filename_in, ...
        'Normalization_energy',Elem.Normalization_energy, ...
        'Nslice',Elem.Nslice, ...
        'Length',Elem.Length, ...
        'xkick',Elem.xkick, ...
        'ykick',Elem.ykick, ...
        'xkick1',Elem.xkick1, ...
        'ykick1',Elem.ykick1, ...
        'xtable',Elem.xtable, ...
        'ytable',Elem.ytable);
    Elem.KickmapStore = struct('default',kickmapData);
    Elem.ActiveKickmap = 'default';
end
