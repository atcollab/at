function varargout = setgirdermisalign(varargin)
% SETGIRDERMISALIGN will add/generate a misalignment structure for THERING that
% will later be used by CALCMISALIGN to generate the actual misalignment
% and apply it to THERING. 
%
% SETGIRDERMISALIGN(FAMNAME1, FAMNAME2, SIGMA, CUTOFF) where FAMNAME1 is
% the name of the element at the start of the girder. It is recommended
% that unique elements be used to identify the girders. SIGMA is a 1x6
% vector specifying [dx, thetax, dy, thetay, ds, thetas] which is the 1
% sigma spread, and CUTOFF an integer specifying how many sigmas before
% excluding the calculated misalignment. If cutoff is 0, then all random
% numbers will be accepted.

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
        indices2 = findcells(THERING, 'FamName',varargin{n});
        if isempty(indices2)
            disp(['Could not find element: ' varargin{n}]);
            return
        else
            famname2 = varargin{n};
        end
    elseif isnumeric(varargin{n})
        % Assumes that if indices are supplied that they are part of the
        % same family of elements.
        indices2 = varargin{n};
        famname2 = THERING{indices2(1)}.FamName;
    else
        error('Family name or index into THERING not in correct format!');
    end
else
    error('Family name or index into THERING not given!');
end

n = n - 1;
if nargin >= n
    if ischar(varargin{n})
        indices1 = findcells(THERING, 'FamName',varargin{n});
        if isempty(indices1)
            disp(['Could not find element: ' varargin{n}]);
            return
        else
            famname1 = varargin{n};
        end
    elseif isnumeric(varargin{n})
        % Assumes that if indices are supplied that they are part of the
        % same family of elements.
        indices1 = varargin{n};
        famname1 = THERING{indices1(1)}.FamName;
    else
        error('Family name or index into THERING not in correct format!');
    end
else
    error('Family name or index into THERING not given!');
end

if length(indices1) ~= length(indices2)
    error('Number of elements specifying the start of the girder does not match those at the end');
else
    for i=1:length(indices1)
        if indices1(i) > indices2(i)
            disp('Seems that the ''start'' of a girder comes AFTER the ''finishing'' element.');
            disp(['   ' famname1 '(' num2str(indices1(i)) ') -> ' famname2 '(' num2str(indices2(i)) ')']);
            return
        else
            girderels{i} = [indices1(i):indices2(i)];
        end
    end
end

mis.ngird = mis.ngird + 1;

mis.gird(mis.ngird).start = famname1;
mis.gird(mis.ngird).finish = famname2;
mis.gird(mis.ngird).ATindex_start = indices1;
mis.gird(mis.ngird).ATindex_finish = indices2;
mis.gird(mis.ngird).ATindex_girder = girderels;
mis.gird(mis.ngird).sigma = sigma;
mis.gird(mis.ngird).cutoff = cutoff;

setappdata(0,'MisalignData',mis);