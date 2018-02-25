
%ATHELP generates the list of Accelerator Toolbox functions

% Generate help file (Contents) use command atupdateContents

ATROOT = getenv('ATROOT');
if isempty(ATROOT)
    ATROOT = atroot;
end

doc(fileparts(ATROOT));
%help(fileparts(ATROOT));
