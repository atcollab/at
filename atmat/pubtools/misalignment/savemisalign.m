function varargout = savemisalign(varargin)
% SAVEMISALIGN will save the current misalignment data structure into a
% file.
%
% SAVEMISALIGN([FILENAME]) where if the optional FILENAME is specified the
% structure will be saved to it else it will open a dialogue box to get the
% user to tell it where and what name to save it to.

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
end

if nargin == 1 & ischar(varargin{1})
    save(varargin{1},'mis');
else
    uisave('mis');
end