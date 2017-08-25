%ATHELP generates the list of Accelerator Toolbox functions
ATROOT = getenv('ATROOT');
if isempty(ATROOT)
    ATROOT = atroot;
end

disp('Physics Tools');
disp('');
help(fullfile(ATROOT,'atphysics'))
disp('Touschek Pivinski');
disp('');
help(fullfile(ATROOT,'atphysics', 'TouschekPiwinski'))
disp('Lattice Tools');
disp('');
help(fullfile(ATROOT, 'lattice'))
help(fullfile(ATROOT, 'atgui'))
disp('Element creation');
disp('');
help(fullfile(ATROOT,'lattice','element_creation'))
disp('atfastring');
disp('');
help(fullfile(ATROOT,'lattice', 'atfastring'))
disp('Survey');
disp('');
help(fullfile(ATROOT,'lattice', 'survey'))
disp('Integrators Tracking Methods');
disp('');
help(fullfile(ATROOT,'..', 'atintegrators'))
help(fullfile(ATROOT,'attrack'))
disp('User defined Integrators');
disp('');
help(fullfile(ATROOT,'..','atintegrators', 'user'))
disp('Matching Tools');
disp('');
help(fullfile(ATROOT,'atmatch'))
disp('Plot Tools to be used with atplot');
disp('');
help(fullfile(ATROOT,'atplot'))
disp('Plot functions to be used with atplot');
help(fullfile(ATROOT,'atplot', 'plotfunctions'))

disp('AT Demos');
disp('');
help(fullfile(ATROOT,'atdemos'))