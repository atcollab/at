function githubsetup(execmode)
%GITHUBSETUP    Setup Matlab for AT tests in GitHib Actions

savepath('pathdef.m');
atmexall;
pyenv('Version',fullfile(getenv('pythonLocation'),'bin','python'),...
      'ExecutioMode', execmode);
end