%ATHELP generates the list of Accelerator Toolbox functions
ATROOT = getenv('ATROOT');
disp('Physics Tools');
disp('');
help([ATROOT,'\atphysics'])
disp('Lattice Tools');
disp('');
help([ATROOT,'\lattice'])
disp('AT Demos');
disp('');
help([ATROOT,'\atdemos'])