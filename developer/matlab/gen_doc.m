function gen_doc()
%GEN_DOC	Build the Matlab documentation
%
%Builds the HTML documentation used by the Matlab help browser and the markdown
%documentation used by Sphinx

w0=warning;  % Save the warning state
warning('off', 'AT:NoRingParam'); % Disable RingParam warnings

[devdir,~, ~]=fileparts(mfilename('fullpath'));
docdir = fullfile(atroot,'..','docs','atdocs','matlab');
sphinxdir = fullfile(atroot,'..','docs','m');
tocfile = fullfile(docdir,'helptoc.xml');
[fid,fmess,fmid] = copyfile(fullfile(devdir,'helptoc.xml'),tocfile);
if fid <= 0
    error(fmid,'%s: %s',tocfile,fmess);
end
fid=openmfile(tocfile,'at');

% Publish top page
src=fullfile(devdir,'mlx','AT_page.mlx');
dst=fullfile(docdir,'AT_page.html');
fprintf('Export %s to %s\n', src, dst);
export(src,dst);

% Getting started
mlxloop('getting_started','Getting started',Run=true)
%sphinxloop('getting_started',Run=true);

% User guide
ugname=fullfile(devdir,'m','ugsummary.m');
hid=openmfile(ugname,'wt');
fprintf(fid,'        <tocitem target="ugsummary.html"\n');
fprintf(fid,'                 image="HelpIcon.USER_GUIDE">\n');
fprintf(fid,'        AT User Guide\n');
fprintf(hid,'%%%% AT User Guide\n%%\n%%%%\n');
%   Loop on UG chapters
for m=atchapters()
    mname = fullfile(devdir,'m', m.id+".m");
    fprintf(hid,'%% <matlab:web(fullfile(docroot,''3ptoolbox'',''atacceleratortoolbox'',''doc'',''%s.html'')) %s>\n%%\n',m.id,m.title);
    fprintf(fid,'            <tocitem target="%s.html">%s</tocitem>\n',m.id,m.title);
    gid=openmfile(mname,'wt');
    fprintf(gid,'%%%% %s\n%% \n%%%%\n', m.title);
    mloop(gid,m.contents);
    fclose(gid);
    publish(mname,'evalCode',false,'outputDir',docdir);
end
fprintf(fid,'        </tocitem>\n');
fclose(hid);
publish(ugname,'evalCode',false,'outputDir',docdir);

% How to...
mlxloop('howtos','How to…');
%sphinxloop('howtos');

% Release notes
mlxloop('release_notes','Release notes');
%sphinxloop('release_notes');

% Web site
fprintf(fid,'        <tocitem target="https://atcollab.github.io/at/" \n');
fprintf(fid,'                 image="$toolbox/matlab/icons/webicon.gif">\n');
fprintf(fid,'        AT Web Site\n');
fprintf(fid,'        </tocitem>\n');

fprintf(fid,'    </tocitem>\n');
fprintf(fid,'</toc>\n');
fclose(fid);

%Build doc search database
try
    builddocsearchdb(docdir);
catch me
    warning( me.message )
end

warning(w0);  % Restore the warning state

    function mloop(fid,mlist)
        for item=mlist
            if startsWith(item,"-")     % Section header
                fprintf(fid,'%%%% %s\n%% \n%%%%\n', eraseBetween(item,1,1));
            elseif startsWith(item,"0") % Plain text
                fprintf(fid,'%% %s\n%%\n',eraseBetween(item,1,1));
            else                        % AT function name
                try
                    h1=h1_line(which(item));
                    did=sprintf('<matlab:doc(''%s'') %s>', h1.name, h1.name);
                    fprintf(fid,'%% %s - %s\n%%\n',did,h1.h1);
                catch err
                    disp(err.message)
                end
            end
        end
    end

    function mlxloop(secdir,secname,varargin)
        sumname=fullfile(devdir,'m',[secdir '.m']);
        sumid=openmfile(sumname,'wt');
        fprintf(fid,'        <tocitem target="%s" \n', [secdir '.html']);
        fprintf(fid,'                 image="$toolbox/matlab/icons/webicon.gif">\n');
        fprintf(fid,'        %s\n',secname);
        fprintf(sumid,'%%%% %s\n%% \n%%%%\n', secname);
        for fn=reshape(dir(fullfile(devdir,secdir,'*.mlx')),1,[])
            source=fullfile(fn.folder,fn.name);
            [~,name,~]=fileparts(fn.name);
            title=regexprep(name,{'^[0-9_ ]*', '_'},{'',' '});
            target=fullfile(secdir,append(name,'.html'));
            targ1=fullfile(docdir,target);
            % Export HTML file for the Matlab help browser
            fprintf('Export %s to %s\n', source, targ1);
            export(source,targ1,varargin{:});
            % Export Markdown file fpr Sphinx
            targ2=fullfile(sphinxdir,secdir,append(name,'.md'));
            fprintf('Export %s to %s\n',source,targ2);
            export(source,targ2,varargin{:});
            fprintf(sumid,'%% <matlab:web(fullfile(docroot,''3ptoolbox'',''atacceleratortoolbox'',''doc'',''%s'')) %s>\n%%\n',target,title);
            fprintf(fid,'            <tocitem target="%s">%s</tocitem>\n',target,title);
        end
        fclose(sumid);
        fprintf('Publish %s to %s\n', sumname, docdir);
        publish(sumname,'evalCode',false,'outputDir',docdir);
        fprintf(fid,'        </tocitem>\n');
        delete(sumname);
    end

    function fid=openmfile(fpath,mode)
        [fid,mess]=fopen(fpath,mode);
        if fid < 0
            error('AT:doc','Cannot open %s: %s',fpath,mess);
        end
    end

end
