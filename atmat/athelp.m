function athelp
%ATHELP generates the list of Accelerator Toolbox functions

% Generate help file (Contents) use command atupdateContents

ATROOT = atroot;

DIR_old = pwd;
cd(fileparts(ATROOT))
doc(cd);
cd(DIR_old)

%help(fileparts(ATROOT));
