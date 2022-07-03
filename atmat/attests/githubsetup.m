function githubsetup(execmode)
%GITHUBSETUP    Setup Matlab for AT tests in GitHib Actions

savepath('pathdef.m');
atmexall;
if ispc
    pyenv("Version",'3.9','ExecutionMode', execmode);
else
    pyenv('Version',fullfile(getenv('pythonLocation'),'bin','python'),...
      'ExecutionMode', execmode);
end
disp(pyenv)
end