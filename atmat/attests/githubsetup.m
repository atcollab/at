function githubsetup(execmode)
%GITHUBSETUP    Setup Matlab for AT tests in GitHib Actions

savepath('pathdef.m');
atmexall;
if ispc
    execfile=fullfile(getenv('pythonLocation'),'pythonw.exe');
else
    execfile=fullfile(getenv('pythonLocation'),'bin','python');
end
pyenv("Version", execfile,'ExecutionMode', execmode);
end