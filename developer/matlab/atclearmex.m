function atclearmex(rootdir)
%ATCLEARMEX Remove all AT mex-files
%   Detailed explanation goes here
% <a href="matlab:web('https://atcollab.github.io/at')">AT Web site</a>.
% For more information, see <a href="matlab: 
% doc('atmat')">atmat</a>.
% For more information, see <a href="matlab: 
% help('atmat')">atmat</a>.

if nargin < 1
    rootdir=atroot;
end
[status,out]=unix(sprintf('find "%s" -name "*.mex*" -exec rm {} \\;', rootdir));
[status,out]=unix(sprintf('rm %s',fullfile(rootdir,'..','atintegrators','*.mex*')));
end