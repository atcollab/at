function githubsetup()
%GITHUBSETUP    Private. Setup Matlab for AT tests in GitHib Actions
%
% This function prepares a workflow in GitHub actions for using AT and
% calling python from Matlab. There is normally no reason to use it in a
% user workflow.

savepath('pathdef.m');
atmexall -fail
if ispc
    execfile=fullfile(getenv('pythonLocation'),'pythonw.exe');
    execmode='InProcess';
elseif ismac
    execfile=fullfile(getenv('pythonLocation'),'bin','python');
    execmode='InProcess';
else
    execfile=fullfile(getenv('pythonLocation'),'bin','python');
    execmode='OutOfProcess';
end
pyenv("Version", execfile,'ExecutionMode', execmode);
disp(pyenv);
end