function atdiag
%ATDIAG - test AT intallation
try
    disp('>> spear3;')
    spear3;
    disp('>> plotbeta');
    plotbeta;
    disp('If you see beta-function plots, AT is installed correctly')
catch
    disp('AT was not installed correctly')
end


