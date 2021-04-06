function atdiag
%ATDIAG Tests AT intallation

try
    disp('>> ring=atradoff(spear3);')
    ring = atradoff(spear3);
    disp('>> atplot(ring);');
    atplot(ring);
    disp('If you see beta-function plots, AT is installed correctly')
catch err
    disp('AT was not installed correctly:');
    disp(err.message);
end


