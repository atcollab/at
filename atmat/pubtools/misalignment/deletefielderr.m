function varargout = deletefielderr(varargin)
% DELETEFIELDERR will delete the field error data structure.

% Note: sould not undo the misalignment here. Deletion will typically be
% used at the end of the study or sometimes at the start to ensure that no
% residual errors are left before commencing the study.

ferr = getappdata(0,'FielderrData');

if isempty(ferr)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
end

rmappdata(0,'FielderrData');

disp('Field error data structure has been deleted.');