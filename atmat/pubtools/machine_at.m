function varargout = machine_at(varargin)

% Machine AT will return the optics of the lattice. Essentially takes what
% twissring returns and restructures it so that its easier to use.
%
% optics = MACHINE_AT([THERING, DP, ELEMENTS, SPECIFIC_OPTICS])
% optics = MACHINE_AT(SPECIFIC_OPTICS)
%
% Defaults to loading THERING from global, DP = 0, uses all elements and
% calculates the dispersion.
% SPECIFIC_OPTICS is a string that defines which particular element to
% return. Eg. if SPECIFIC_OPTICS = 'betax', then MACHINE_AT will only
% return those numbers.
%
% Other options are 'file' which will allow one to export the optical
% functions to an xls spreadsheet and 'line' to specify if the THERING used
% is actually a transfer line. If using the transfer line then the user
% must also provide the initial conditions as specified by TWISSLINE.
%
% Added export of element names as well 14/11/06
% Added export to xls option. 12/05/05
% Eugene 28/05/04

nspec_opt = 0;
exportfile = 0;
transferline = 0;
tdin = [];
for i=1:nargin
    if ischar(varargin{i}) && strcmpi(varargin{i},'file')
        exportfile = 1;
    elseif ischar(varargin{i}) && strcmpi(varargin{i},'line')
        transferline = 1;
    elseif isstruct(varargin{i})
        % Of all inputs there is only one struct type and that is for
        % twissline.
        tdin = varargin{i};
    else
        if isstr(varargin{i})
            nspec_opt = nspec_opt + 1;
            spec_opt{nspec_opt} = varargin{i};
        end
    end
end
        
args = 0;
args = args + 1;
if nargin >= args && ~isstr(varargin{args}) && ~isstruct(varargin{args})
    line_ring = varargin{args};
else
    global THERING
    line_ring = THERING;
end

args = args + 1;
if nargin >= args && ~isstr(varargin{args}) && ~isstruct(varargin{args})
    dp = varargin{args};
else
    dp = 0;
end

args = args + 1;
if nargin >= args && ~isstr(varargin{args}) && ~isstruct(varargin{args})
    elements = varargin{args};
else
    elements = 1:length(line_ring)+1;
end

% Check that input struct supplied of 'line' is used
if transferline
    if isempty(tdin)
        error('User must provide the twiss data input structure');
    end
    TD = twissline(line_ring, dp, tdin, elements, 'chrom', 1e-6);
else
    TD = twissring(line_ring, dp, elements, 'chrom', 1e-6);
end

% Group element names into a cell array
for i=1:length(elements)
    % Circular indexing
    iind = mod(elements(i)-1,length(line_ring))+1;
    elemnames{i,1} = line_ring{iind}.FamName;
    if isfield(line_ring{iind},'Length')
        elemLeff(i,1) = line_ring{iind}.Length;
    else
        elemLeff(i,1) = 0;
    end
end
optics.elemnames = elemnames;
optics.elemLeff = elemLeff;

temp = cat(1, TD.beta);
optics.betax = temp(:,1);
optics.betay = temp(:,2);

temp = cat(1, TD.alpha);
optics.alphax = temp(:,1);
optics.alphay = temp(:,2);

temp = cat(2, TD.Dispersion);
optics.etax = temp(1,:)';
optics.etapx = temp(2,:)';
optics.etay = temp(3,:)';
optics.etapy = temp(4,:)';

temp = cat(2, TD.ClosedOrbit);
optics.x = temp(1,:)';
optics.px = temp(2,:)';
optics.y = temp(3,:)';
optics.py = temp(4,:)';

temp = cat(1,TD.mu);
optics.nux = temp(:,1)/(2*pi);
optics.nuy = temp(:,2)/(2*pi);

optics.spos = cat(1,TD.SPos);

if nspec_opt >= 1
    varargout{1} = optics.(spec_opt{1});
else
    varargout{1} = optics;
end

if exportfile
    [filename pathname] = uiputfile('*.xls','Excel spreadsheet');
    entrystr = fieldnames(optics);
    temp = optics.(entrystr{1});
    for i=2:length(entrystr)
        temp = cat(2,temp,optics.(entrystr{i}));
    end
    xlswrite([pathname filename],entrystr','opticalparam','A1');
    xlswrite([pathname filename],temp,'opticalparam','A2');
end
