function atm2html(varargin)
%MAKEDOC_HTML - Generate new MML, SOLEIL and AT HTML help files
%  makedoc_html
%
%  HOWTO
%  1. Make sure to update and run toolboxUpdateHeader.m
%  2. Update history.txt appropriately, including w current version
%  3. Update overview.html file with the version/date/link to zip:
%     edit external/m2html/templates/at/about.html
%  4. Need to install graphviz fro graph dependency
%     see: https://graphviz.org/

%
%%  Written  by Laurent S. Nadolski

Options = {...
    'htmldir','doc_html', ...
    'indexfile', 'atindex',...
    'recursive','on', ...
    'graph','on', ...
    'todo', 'on', ...
    'globalHypertextLinks', 'on', ...
    'helptocxml', 'on', ...
    'template', 'frame'... % template for AT
    'index','menu', ... % this part in mandatory with frame
    'global', 'on', ...
    'save', 'on', ...
    'download','on', ...
    'search', 'off', ... % does not work properly
    'verbose', 'on', ...
    'syntaxHighlighting', 'on', ...
    };
%    'ignoredDir', {{'.svn' 'cvs'}}, ...

DirectoryStart = pwd;

Directory = atroot;

cd(Directory);

% Make AT HTML help
% Remove old doc
cd ..

try
    if isfolder('doc_html')
        rmdir('doc_html','s');
    end
catch
    disp('rmdir error')
end

DirectoryList = {...
    'atmat', ...
    'atintegrators',...
    'machine_data',...
    };

m2html('mfiles', DirectoryList, Options{:});

% Move doc_html directory to AT
%movefile('doc_html', 'at');

cd(DirectoryStart);


end
