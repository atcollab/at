function varargout = setmisalign(varargin)
% SETMISALIGN will add/generate a misalignment structure for THERING that
% will later be used by CALCMISALIGN to generate the actual misalignment
% and apply it to THERING. 
%
% SETMISALIGN(FAMNAME, SIGMA, CUTOFF[, LATTICE]) where FAMNAME is the name 
% of the family of elements, SIGMA is a 1x6 vector specifying [dx, thetax, 
% dy, thetay, ds, thetas] which is the 1 sigma spread, and CUTOFF an integer
% specifying how many sigmas before excluding the calculated misalignment.
% If cutoff is 0, then all random numbers will be accepted.

global THERING

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('Generating the first misalignemnt datastructure now');
    mis.date = datestr(now,31);
    mis.nind = 0;
    mis.ngird = 0;
end

% Parse input data
n = 4;
if nargin >= n & iscell(varargin{n})
    LATTICE = varargin{n};
else
    disp('Using THERING');
    LATTICE = THERING;
end

n = n - 1;
if nargin >= n & isnumeric(varargin{n})
    cutoff = varargin{n};
else
    disp('Defaulting cutoff to 3.');
    cutoff = 3;
end

n = n - 1;
if nargin >= n & isnumeric(varargin{n}) & length(varargin{n}) == 6
    sigma = varargin{n};
else
    error('One sigma spread of the misalignment not specified properly!');
end

n = n - 1;
if nargin >= n
    if ischar(varargin{n})
        indices = findcells(LATTICE, 'FamName',varargin{n});
        famname = varargin{n};
    elseif isnumeric(varargin{n})
        % Assumes that if indices are supplied that they are part of the
        % same family of elements.
        indices = varargin{n};
        famname = LATTICE{indices(1)}.FamName;
    elseif iscell(varargin{n})
        for i=1:length(varargin{n})
            setmisalign(varargin{n}{i},sigma,cutoff);
        end
        return
    else
        error('Family name or index into THERING not in correct format!');
    end
else
    error('Family name or index into THERING not given!');
end

mis.nind = mis.nind + 1;

mis.ind(mis.nind).name = famname;
mis.ind(mis.nind).ATindex = indices;
mis.ind(mis.nind).sigma = sigma;
mis.ind(mis.nind).cutoff = cutoff;

setappdata(0,'MisalignData',mis);