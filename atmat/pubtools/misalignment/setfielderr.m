function varargout = setfielderr(varargin)
% SETFIELDERR will add/generate a field error structure for THERING that
% will later be used by CALCFIELDERR to generate the actual field error
% and apply it to THERING. 
%
% SETFIELDERR(FAMNAME, SIGMA_A, SIGMA_B, CUTOFF) where FAMNAME is the name
% of the family of elements, SIGMA_A/SIGMA_B is a scalar or vector
% specifying 1 sigma spread of the field or multipole fields, polynomA and
% polynomB respectively. CUTOFF an integer specifying how many sigmas
% before excluding the calculated misalignment. If cutoff is 0, then all
% random numbers will be accepted.

global THERING

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

ferr = getappdata(0,'FielderrData');

if isempty(ferr)
    disp('Generating the first field error datastructure now');
    ferr.date = datestr(now,31);
    ferr.CreatedBy = 'setfielderr';
    ferr.nind = 0;
    ferr.ngird = 0;
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
if nargin >= n & isnumeric(varargin{n})
    sigmaB = varargin{n};
else
    error('One sigma spread of the field error not specified properly! (SIGMA_B)');
end

n = n - 1;
if nargin >= n & isnumeric(varargin{n})
    sigmaA = varargin{n};
else
    error('One sigma spread of the field error not specified properly! (SIGMA_A)');
end

n = n - 1;
if nargin >= n
    if ischar(varargin{n})
        indices = findcells(THERING, 'FamName',varargin{n});
        famname = varargin{n};
    elseif isnumeric(varargin{n})
        % Assumes that if indices are supplied that they are part of the
        % same family of elements.
        indices = varargin{n};
        famname = THERING{indices(1)}.FamName;
    elseif iscell(varargin{n})
        for i=1:length(varargin{n})
            setfielderr(varargin{n}{i},sigma,cutoff);
        end
        return
    else
        error('Family name or index into THERING not in correct format!');
    end
else
    error('Family name or index into THERING not given!');
end

ferr.nind = ferr.nind + 1;

ferr.ind(ferr.nind).name = famname;
ferr.ind(ferr.nind).ATindex = indices;
ferr.ind(ferr.nind).sigmaA = sigmaA;
ferr.ind(ferr.nind).sigmaB = sigmaB;
ferr.ind(ferr.nind).cutoff = cutoff;

setappdata(0,'FielderrData',ferr);