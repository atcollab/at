function githubsetup()
%GITHUBSETUP    Setup Matlab for AT tests in GitHib Actions

savepath('pathdef.m');
atmexall;
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