function varargout = loadmisalign(varargin)
% LOADMISALIGN will load the misalignment data structure from a file that
% was previously saved.
%
% LOADMISALIGN([FILENAME]) where if the FILENAME is specified, it will try
% to load the file. Otherwise it will open a dialogue for the user to pick
% the file to load.

if nargin == 1 & ischar(varargin{1})
    load(varargin{1});
else
    uiload;
end

if exist('mis','var') & isstruct(mis) & isfield(mis,'nind')
    setappdata(0,'MisalignData',mis);
else
    error('File being loaded does not seem to contain the misalignment data structure');
end