%ATHELP generates the list of Accelerator Toolbox functions
ATROOT = getenv('ATROOT');

disp('Physics Tools');
disp('');
help(fullfile(ATROOT,'atphysics'))
disp('Lattice Tools');
disp('');
help(fullfile(ATROOT, 'lattice'))
disp('AT Demos');
disp('');
help(fullfile(ATROOT,'atdemos'))