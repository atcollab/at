function varargout = applymisalign(varargin)
% APPLYMISALIGN will read the misalignment data structure and apply it to
% THERING elements. 
%
% LATTICE = APPLYMISALIGN([SEED, LATTICE]). If more than one seed has been 
% generated using CALCMISALIGN then repeated calling of APPLYMISALIGN will 
% apply the subsequent seeds untill all have been used and the user will 
% have to generate more seeds. The seed that is currently being used is 
% saved in the misalignment datastructure under the field 'currseed'. If 
% SEED is defined then that particular seed will be applied and the running
% 'currseed' will not be incremented.
%

global THERING

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

[reg prop] = parseparams(varargin);
verbatim = 1;
for i=1:length(prop)
    if strcmpi(prop,'quiet')
        verbatim = 0;
    end
end

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
end

n = 2;
if nargin >= n & iscell(varargin{n})
    LATTICE = varargin{n};
else
    disp('Using THERING');
    LATTICE = THERING;
end

n = n-1;
if nargin >= n & isnumeric(varargin{n})
    mis.currseed = varargin{n};
else
    mis.currseed = mis.currseed + 1;
    if mis.currseed > mis.numseeds
        disp('All seeds used. Doing nothing')
        return
    end
end
    
% Indices of elements where misalignemnt data may need to be applied.
indices = find(cell2mat({mis.data.used}));
% Cancatenate all the misalignment data into one huge matrix.
temp = cell2mat({mis.data(indices).val});
% Filtering out the seeds being used.
misalignval = zeros(6,length(indices));
misalignval = temp(1:6,[mis.currseed:mis.numseeds:end]);

% Separate out the individual shifts
xshift = misalignval(1,:);
xrot   = misalignval(2,:);
yshift = misalignval(3,:);
yrot   = misalignval(4,:);
sshift = misalignval(5,:);
srot   = misalignval(6,:);

% Hack so that I don't have to change the functions below to accept
% arbritrary lattices
TEMP = THERING;
THERING = LATTICE;

% Set and add the rotations because these affect T1 and T2.
setshift(indices,xshift,yshift);
addxrot(indices,xrot);
addyrot(indices,yrot);
% This is ok to set since the roll affects R1 and R2 only.
addsrot(indices,srot);
% INCOMPLETE
% addsshift(indices,sshift);

% Put everything back
LATTICE = THERING;
THERING = TEMP;
clear TEMP;

if verbatim
    disp(['Applied random misalignemnt. Using seed number ' num2str(mis.currseed)]);
end

% Saving some info for plotting and analysis of applied misalignments
% get the SPos data
TD3 = twissring(LATTICE,0,1:length(LATTICE));
S3 = cat(2,TD3.SPos);
mis.applied(mis.currseed).spos = S3;
mis.applied(mis.currseed).deltax = zeros(1,length(LATTICE)); 
mis.applied(mis.currseed).deltax(indices) = xshift;
mis.applied(mis.currseed).deltay = zeros(1,length(LATTICE)); 
mis.applied(mis.currseed).deltay(indices) = yshift;
mis.applied(mis.currseed).deltas = zeros(1,length(LATTICE)); 
mis.applied(mis.currseed).deltas(indices) = sshift;


setappdata(0,'MisalignData',mis);

if nargout == 1
    varargout{1} = LATTICE;
end