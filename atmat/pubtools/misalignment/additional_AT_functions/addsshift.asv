function addsshift(elemindex, ds)
% ADDSSHIFT will translate the element along the 's' vector by extending
% and contracting drifts that exist on either side of the element. If there
% are other nondrift elements beside the elements (eg split elements) there
% will be an error and the user will be prompted to check the element they
% wish to move and ensure that it has adjacent drifts (ie space to move).
% The program will prescan "elemindex" and will group elements if the
% indices are sequential. Example: 1 3 5 6 8 10 -> will be grouped such
% that 5 and 6 are considered a single element.
%
% !!! WARNING !!! In order to undo the changes one has to reload to lattice
% file as the shifts are not stored here and will not know what to do to
% undo the length changes. If using the misalignment tools then the
% UNDOMISALIGN function should be able to undo the changes that were
% applied to it.
%
% See also ADDXROT, ADDYROT, ADDSROT, ADDSHIFT, SETXROT, SETYROT, SETSROT

global THERING FAMLIST

% Elemindex and ds have to be the same size
if length(elemindex) ~= length(ds)
    disp('Number of elements and number of shifts are not the same');
    disp('Nothing done');
    return
end

% Find and group sequential elements
ind = [];
if length(elemindex) > 1
    temp = elemindex(2:end) - elemindex(1:end-1);
    ind = find(temp < 2);
end

% Check that sequential elements are being shifted by the same amount.
for i=1:length(ind)
    if ds(ind(i)) ~= ds(ind(i)+1)
        disp('Two sequential elements (ie no drift between them) are moved by different amounts.')
        disp('!! This cannot be done !!');
        disp([' -> take a look at element index: ' num2str(ind(i)) ' and ' num2str(ind(i) + 1)]);
        return
    end
end

% Index the "front" and "back" element and check that they are drifts




return
frontdrift = 0;
backdrift = 0;
for i=1:length(elemindex)
    currel = elemindex(i);

    % Find first drift in front of element
    ind = currel;
    while frontdrift == 0
        if isempty(regexpi(THERING{ind}.PassMethod,'drift')) | THERING{ind}.Length == 0
            % Assume circular ring
            if ind == length(THERING)
                ind = 1;
            else
                ind = ind + 1;
            end
        else
            frontdrift = ind;
        end
    end

    % Find first drift behind the element
    ind = currel;
    while backdrift == 0
        if isempty(regexpi(THERING{ind}.PassMethod,'drift')) | THERING{ind}.Length == 0
            % Assume circular ring
            if ind == 1
                ind = length(THERING);
            else
                ind = ind - 1;
            end
        else
            backdrift = ind;
        end
    end

    THERING{backdrift}.Length = THERING{backdrift}.Length + ds(i);
    THERING{frontdrift}.Length = THERING{frontdrift}.Length + ds(i);
end