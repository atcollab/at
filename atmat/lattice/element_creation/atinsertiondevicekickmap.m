function Elem = atinsertiondevicekickmap(fname,length,normalizationEnergy,varargin)
%ATINSERTIONDEVICEKICKMAP Create an insertion device kick-map element
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME)
%   Create an empty zero-length insertion device normalized at 0 GeV, with
%   PASSMETHOD='DriftPass'.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,LENGTH,NORMALIZATION_ENERGY,...
%   [PASSMETHOD],'FIELD1',VALUE1,...)
%   Create an element from field/value pairs. PASSMETHOD defaults to
%   'DriftPass'. FAMNAME, LENGTH and NORMALIZATION_ENERGY are positional;
%   all remaining element attributes are supplied as field/value pairs.
%
% Inputs:
%   FAMNAME              Family name
%   LENGTH               Insertion device length [m]
%   NORMALIZATION_ENERGY Energy in GeV used to normalize the kick map
%   PASSMETHOD           Tracking function. Default: 'DriftPass'
%
% Field/value pairs:
%   'Filename_in'        Source kick-map file name
%   'Nslice'             Number of integration slices
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
%   id = atinsertiondevicekickmap('ID',2.0,2.75,'IdTablePass',...
%       'Filename_in','id_kicks.txt','Nslice',25,...
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

if nargin < 2 || isempty(length)
    length = 0.0;
end
if nargin < 3 || isempty(normalizationEnergy)
    normalizationEnergy = 0.0;
end

[rsrc,method] = decodeatargs({'DriftPass'},varargin);
Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap', ...
    'Length',length,'Normalization_energy',normalizationEnergy,rsrc{:});
