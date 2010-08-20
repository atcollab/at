function varargout = applyfielderr(varargin)
% APPLYFIELDERR will read the field error data structure and apply it to
% THERING elements. 
%
% APPLYFIELDERR([SEED],['quiet']). If more than one seed has been generated using
% CALCFIELDERR then repeated calling of APPLYFIELDERR will apply the
% subsequent seeds untill all have been used and the user will have to
% generate more seeds. The seed that is currently being used is saved in
% the field error datastructure under the field 'currseed'. If SEED is
% defined then that particular seed will be applied and the running
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

ferr = getappdata(0,'FielderrData');

if isempty(ferr)
    disp('No field error data found. See SETFIELDERR for more info');
    return
end

if nargin == 1 & isnumeric(varargin{1})
    ferr.currseed = varargin{1};
else
    ferr.currseed = ferr.currseed + 1;
    if ferr.currseed > ferr.numseeds
        disp('All seeds used. Doing nothing')
        return
    end
end
    
% Indices of elements where field error data may need to be applied.
indices = find(cell2mat({ferr.data.used}));
for i=indices
    THERING{i}.PolynomA = THERING{i}.PolynomA + ferr.data(i).valA(:,ferr.currseed)';
    THERING{i}.PolynomB = THERING{i}.PolynomB + ferr.data(i).valB(:,ferr.currseed)';
    if isfield(THERING{i},'K')
        THERING{i}.K = THERING{i}.PolynomB(2);
    end
end

if verbatim
    disp(['Applied random field error. Using seed number ' num2str(ferr.currseed)]);
end

% Saving some info for plotting and analysis of applied field errors
% get the SPos data
% TD3 = twissring(THERING,0,1:length(THERING));
% S3 = cat(2,TD3.SPos);
% ferr.applied(ferr.currseed).spos = S3;
% ferr.applied(ferr.currseed).deltax = zeros(1,length(THERING)); 
% ferr.applied(ferr.currseed).deltax(indices) = xshift;
% ferr.applied(ferr.currseed).deltay = zeros(1,length(THERING)); 
% ferr.applied(ferr.currseed).deltay(indices) = yshift;
% ferr.applied(ferr.currseed).deltas = zeros(1,length(THERING)); 
% ferr.applied(ferr.currseed).deltas(indices) = sshift;


setappdata(0,'FielderrData',ferr);
