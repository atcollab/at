function Elem = atinsertiondevicekickmap(fname,varargin)
% atinsertiondevicekickmap creates an insertion device kick-map element
% Elem = atinsertiondevicekickmap(FAMNAME,[PASSMETHOD],[KEY,VALUE]...)
%
% FAMNAME       Family name
% PASSMETHOD    Tracking function. Default: 'DriftPass'
% KEY, VALUE    Element attributes
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

hasmethod = ~isempty(varargin) && ...
    (ischar(varargin{1}) || isstring(varargin{1})) && ...
    endsWith(varargin{1},'Pass');
optionstart = 1 + hasmethod;
options = varargin(optionstart:end);
namevalue = mod(length(options),2) == 0 && ...
    all(cellfun(@(arg) ischar(arg) || isstring(arg),options(1:2:end)));
legacy = length(varargin) >= 11 && ~namevalue;

if legacy
    method = varargin{1};
    fields = {'Filename_in',varargin{2}, ...
        'Normalization_energy',varargin{3},'Nslice',varargin{4}, ...
        'Length',varargin{5},'xkick',varargin{6},'ykick',varargin{7}, ...
        'xkick1',varargin{8},'ykick1',varargin{9}, ...
        'xtable',varargin{10},'ytable',varargin{11}};
    Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap', ...
        fields{:},varargin{12:end});
else
    [rsrc,method] = decodeatargs({'DriftPass'},varargin);
    Elem = atbaselem(fname,method,'Class','InsertionDeviceKickMap',rsrc{:});
end
