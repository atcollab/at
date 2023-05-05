function setversion(version,release)
%SETVERSION Set the version of AT
%
%SETVERSION()
%   Keep the version number and try to get a release string from the "pull
%   request" number. If not available, the release is set to 'atcollab'.
%
%SETVERSION(VERSION)
%   Set a new version number. The release is set to 'atcollab'.
%
%SETVERSION(VERSION,RELEASE)
%   Set the version and release.

% Try to get the Pull Request number
if nargin < 2
    release='atcollab';
end
if nargin < 1
    [err,outp]=unix('/usr/local/bin/gh pr status --json number --jq .currentBranch.number');
    if err == 0
        prnumber=str2double(outp(1:end-1));  % skip ending newline
        if isfinite(prnumber)
            release=sprintf('#%i',prnumber);
        else
            error('No pull request corresponding to the current branch');
        end
    else
        error('%s\n%s','The release cannot be obtained from the pull request number.', ...
            'Check the "gh" command');
    end
end

versfile=fullfile(atroot,'Contents.m');
tmpfile=fullfile(tempdir,'Contents.m');
fout=fopen(tmpfile,'wt');
fin=fopen(versfile,'rt');
line=fgetl(fin); %#ok<NASGU>                    Header: accelerator toolbox

vv=textscan(fgetl(fin), '%*s %*s %s %s %s');  % Version string
if nargin < 1 || isempty(version)
    version=vv{1}{1};
end
fprintf(fout,'%% Accelerator Toolbox\n');
fprintf(fout,'%% Version %s (%s) %s\n',version,release,datetime('today'));

line=fgetl(fin);
while ~isnumeric(line)
    fprintf(fout,'%s\n',line);
    line=fgetl(fin);
end
fclose(fin);
fclose(fout);
[success,message,messageid]=copyfile(tmpfile,versfile);
if ~success
    error(messageid,message);
end
delete(tmpfile);
gen_help();
end