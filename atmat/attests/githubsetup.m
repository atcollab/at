function githubsetup(varargin)
%GITHUBSETUP    Private. Setup Matlab for AT tests in GitHib Actions
%
% This function prepares a workflow in GitHub actions for using AT and
% calling python from Matlab. There is normally no reason to use it in a
% user workflow.

savepath('pathdef.m');
mex -setup
mex -setup C++
atmexall('-fail', varargin{:});
if ispc
    execfile=fullfile(getenv('pythonLocation'),'pythonw.exe');
else
    execfile=fullfile(getenv('pythonLocation'),'bin','python');
end
pyenv("Version", execfile,'ExecutionMode', 'InProcess');
disp(pyenv);
end