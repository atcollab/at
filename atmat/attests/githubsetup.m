function githubsetup(execmode)
%GITHUBSETUP    Setup Matlab for AT tests in GitHib Actions

savepath('pathdef.m');
atmexall -c_only;
pyenv('Version',fullfile(getenv('pythonLocation'),'bin','python'),...
      'ExecutionMode', execmode);
end