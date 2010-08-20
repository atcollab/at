function varargout = deletemisalign(varargin)
% DELETEMISALIGN will delete the misalignment data structure.

% Note: sould not undo the misalignment here. Deletion will typically be
% used at the end of the study or sometimes at the start to ensure that no
% residual errors are left before commencing the study.

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
end

% try
%     undomisalign;
% catch
%     disp('WARNING! Could not undo misalign. THERING may have been reloaded')
%     disp(' or another problem has occured. Reload lattice before continuing');
% end
rmappdata(0,'MisalignData');

disp('Misalignment data structure has been deleted.');