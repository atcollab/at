function gen_toc()
%GEN_TOC	Build the HTML files used by the Matlab help browser

[here,~, ~]=fileparts(mfilename('fullpath'));
docdir = fullfile(atroot,'..','docs','matlab');
tocfile = fullfile(docdir,'helptoc.xml');
[fid,fmess,fmid] = copyfile(fullfile(here,'helptoc.xml'),tocfile);
if fid <= 0
    error(fmid,'%s: %s',tocfile,fmess);
end
ugname=fullfile('m','ugsummary.m');
fid=openmfile(tocfile,'at');
hid=openmfile(ugname,'wt');

fprintf(fid,'        <tocitem target="ugsummary.html "\n');
fprintf(fid,'            image="HelpIcon.USER_GUIDE">AT User Guide\n');
fprintf(hid,'%%%% AT User Guide\n%%\n%%%%\n');
for m=atchapters()
    mname = fullfile('m', m.id+".m");
    fprintf(hid,'%% <matlab:web(fullfile(docroot,''3ptoolbox'',''atacceleratortoolbox'',''doc'',''%s.html'')) %s>\n%%\n',m.id,m.title);
    fprintf(fid,'            <tocitem target="%s.html">%s</tocitem>\n',m.id,m.title);
    gid=openmfile(mname,'wt');
    fprintf(gid,'%%%% %s\n%% \n%%%%\n', m.title);
    mloop(gid,m.contents);
    fclose(gid);
    publish(mname,'evalCode',false,'outputDir',docdir);
end
fprintf(fid,'        </tocitem>\n');

fprintf(fid,'        <tocitem target="https://atcollab.github.io/at/" \n');
fprintf(fid,'                 image="$toolbox/matlab/icons/webicon.gif">\n');
fprintf(fid,'        AT Web Site\n');
fprintf(fid,'        </tocitem>\n');
fprintf(fid,'    </tocitem>\n');
fprintf(fid,'</toc>\n');
fclose(fid);
fclose(hid);
publish(ugname,'evalCode',false,'outputDir',docdir);

    function mloop(fid,mlist)
        for item=mlist
            if startsWith(item,"-")     % Section header
                fprintf(fid,'%%%% %s\n%% \n%%%%\n', eraseBetween(item,1,1));
            elseif startsWith(item,"0") % Plain text
                fprintf(fid,'%% %s\n',eraseBetween(item,1,1));
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

    function fid=openmfile(fpath,mode)
        [fid,mess]=fopen(fpath,mode);
        if fid < 0
            error('AT:doc','Cannot open %s: %s',fpath,mess);
        end
    end

end
