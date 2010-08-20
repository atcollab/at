function varargout = undomisalign(varargin)
% UNDOMISALIGN will 'zero' all misalignments in THERING 

global THERING

if ~exist('THERING','var')
    disp('Please load the model first');
    return
end

mis = getappdata(0,'MisalignData');

if isempty(mis)
    disp('No misalignment data found. See SETMISALIGN for more info');
    return
end


if isfield(mis,'data') & mis.currseed > 0
    % Indices of elements where misalignemnt data may need to be applied.
    indices = find(cell2mat({mis.data.used}));
    % Cancatenate all the misalignment data into one huge matrix.
    temp = cell2mat({mis.data(indices).val});
    % Filtering out the seeds being used.
    misalignval = zeros(6,length(indices));
    misalignval = temp(1:6,[mis.currseed:mis.numseeds:end]);

    % Separate out the individual shifts.
    % !! SSHIFT is unique in that to remove/reset the individual shifts
    % must be undone and not just set to zero as this involves actual
    % changes to the drifts.
    xshift = -misalignval(1,:);
    xrot   = -misalignval(2,:);
    yshift = -misalignval(3,:);
    yrot   = -misalignval(4,:);
    sshift = -misalignval(5,:);
    srot   = -misalignval(6,:);

    % Set and add the rotations because these affect T1 and T2.
    addshift(indices,xshift,yshift);
    addxrot(indices,xrot);
    addyrot(indices,yrot);
    addsrot(indices,srot);
    % addsshift(indices,sshift);
end

disp('Misalignments have been reset to zero');

