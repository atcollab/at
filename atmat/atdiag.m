function atdiag
%ATDIAG Tests AT intallation

try
    disp('>> ring=atradoff(spear3);')
    ring = atradoff(spear3);
    disp('>> atplot(ring);');
    atplot(ring);
    disp('If you see beta-function plots, AT is installed correctly')
catch err
    if strcmp(err.identifier,'at:missingMex')
        fprintf('AT is not ready.\nCall "atmexall()" to compile the sources.\n')
    else
        fprintf('AT is not installed correctly:\n%s.\n', err.message)
    end
end


