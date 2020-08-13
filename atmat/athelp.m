function athelp(varargin)
%ATHELP generate the list of Accelerator Toolbox functions
%
%  INPUTS
%    No argument - open the help file in Matlab browser
%    'new'       - force the update of the documentation, which requires a few more seconds
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

NEWFLAG=getflag(varargin,'new');

ATROOT = atroot;

fileName = 'athelp.mat';

DIR_old = pwd;
cd(fileparts(ATROOT))

%for compatibility with previous Matlab version 2013b for example
%folder should not include full path
[~, folder] = fileparts(cd);

if ~isfile(fileName) || NEWFLAG
    fprintf('** Generating doc for a few seconds **\n');
    html = help2html(folder,'AT DOC','-doc');
    save(fileName,'html')
else
    d = dir(fileName);
    fprintf('Last version of the file %s \n To update the doc issue athelp(''new'')\n', d.date);
    load(fileName, 'html');
end

% Display documentation
web(['text://' html], '-helpbrowser');

cd(DIR_old)
