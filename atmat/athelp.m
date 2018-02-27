function athelp
%ATHELP generates the list of Accelerator Toolbox functions
%
%  EXAMPLES
%  1. athelp: full help.
%  2. for selected help, use help directory where directory is
%      help atintegrators
%      help atmat
%      help atdemos
%      help atgui
%      help atmatch
%      help atphysics
%      help linearoptics
%      help longitudinaldynamics
%      help nonlineardynamics
%      help atplot
%      help plotfunctions
%      help lattice
%      help element_creation
%      help pubtools
%      help survey
%      help lattice_tools
%      help LatticeTuningFunctions
%      help machine_date
%      help tuneandchromaticity
%      help touschekpiwinski
%      help radiation
%      help parametersummaryfunctions
%           
%
%  See also help

% Generate help file (Contents) use command atupdateContents

ATROOT = atroot;

DIR_old = pwd;
cd(fileparts(ATROOT))

%for comaptbility with previous Matlab version 2013b for example
%folder should not include full path
[~, folder] = fileparts(cd);
doc(folder);

cd(DIR_old)

%help(fileparts(ATROOT));
