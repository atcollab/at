function atdiag
%ATDIAG Tests AT intallation

try
    disp('>> spear3;')
    ring = spear3;
    disp('>> plotbeta');
    atplot(ring);
    disp('If you see beta-function plots, AT is installed correctly')
catch
    disp('AT was not installed correctly')
end


