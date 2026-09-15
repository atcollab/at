function Elem = atinsertiondevicekickmap(fname,varargin)
%ATINSERTIONDEVICEKICKMAP Create an insertion device kick-map element
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME)
%   Create an empty insertion device with PASSMETHOD='DriftPass'.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,PASSMETHOD,'FIELD1',VALUE1,...)
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,'FIELD1',VALUE1,...)
%   Create an element from field/value pairs. PASSMETHOD defaults to
%   'DriftPass'. This is the standard AT syntax and is used for elements
%   written to an M-file by pyAT. Any element field may be supplied,
%   including KickmapStore and ActiveKickmap.
%
% ELEM=ATINSERTIONDEVICEKICKMAP(FAMNAME,PASSMETHOD,FILENAME,...
%   NORMALIZATION_ENERGY,NSLICE,LENGTH,XKICK,YKICK,XKICK1,YKICK1,...
%   XTABLE,YTABLE)
%   Create an element with the original positional syntax. This form is
%   retained for compatibility with existing MATLAB code.
%
% Inputs:
%   FAMNAME              Family name
%   PASSMETHOD           Tracking function. Default: 'DriftPass'
%   FILENAME             Source kick-map file name
%   NORMALIZATION_ENERGY Energy in GeV used to normalize the kick map
%   NSLICE               Number of integration slices
%   LENGTH               Insertion device length [m]
%   XKICK, YKICK         Second-order horizontal and vertical kick tables
%   XKICK1, YKICK1       First-order horizontal and vertical kick tables
%   XTABLE, YTABLE       Horizontal and vertical table coordinates
%
% Examples:
%   emptyid = atinsertiondevicekickmap('ID');
%   id = atinsertiondevicekickmap('ID','IdTablePass',...
%       'Filename_in','id_kicks.txt','Normalization_energy',2.75,...
%       'Nslice',25,'Length',2.0,'xkick',xkick,'ykick',ykick,...
%       'xkick1',xkick1,'ykick1',ykick1,'xtable',x,'ytable',y);
%   id = atinsertiondevicekickmap('ID','IdTablePass','id_kicks.txt',...
%       2.75,25,2.0,xkick,ykick,xkick1,ykick1,x,y);
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

% The original constructor used fixed positional arguments. New M-files use
% the usual AT field/value syntax so that additional fields, such as the
% multi-kickmap store, can be preserved without changing this signature.
hasPositionalPassMethod = ~isempty(varargin) && ...
    (ischar(varargin{1}) || isstring(varargin{1})) && ...
    endsWith(varargin{1},'Pass');
firstAttribute = 1 + hasPositionalPassMethod;
attributeArguments = varargin(firstAttribute:end);
usesNameValueSyntax = mod(length(attributeArguments),2) == 0 && ...
    all(cellfun(@(arg) ischar(arg) || isstring(arg), ...
    attributeArguments(1:2:end)));
usesPositionalCompatibilitySyntax = hasPositionalPassMethod && ...
    length(varargin) >= 11 && ~usesNameValueSyntax;

if usesPositionalCompatibilitySyntax
    % Map the original fixed argument order to named element fields.
    method = varargin{1};
    fields = {'Filename_in',varargin{2}, ...
        'Normalization_energy',varargin{3},'Nslice',varargin{4}, ...
        'Length',varargin{5},'xkick',varargin{6},'ykick',varargin{7}, ...
        'xkick1',varargin{8},'ykick1',varargin{9}, ...
        'xtable',varargin{10},'ytable',varargin{11}};
    Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap', ...
        fields{:},varargin{12:end});
else
    % Standard AT constructor: optional pass method followed by field/value
    % pairs. atbaselem also accepts PassMethod and FamName overrides.
    [rsrc,method] = decodeatargs({'DriftPass'},varargin);
    Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap',rsrc{:});
end
