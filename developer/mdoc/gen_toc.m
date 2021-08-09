function gen_toc()
%GEN_TOC	Build the HTML files used by the Matlab help browser

[here,~, ~]=fileparts(mfilename('fullpath'));
docdir = fullfile(atroot,'..','doc_matlab');
tocfile = fullfile(docdir,'helptoc.xml');
[fid,fmess,fmid] = copyfile(fullfile(here,'helptoc.xml'),tocfile);
if fid <= 0
    error(fmid,'%s: %s',tocfile,fmess);
end
[fid,fmess]=fopen(tocfile,'at');
if fid < 0
    error('AT:doc','%s: %s',tocfile,fmess);
end

fprintf(fid,'        <tocitem target="mytbx_ug_intro.html "\n');
fprintf(fid,'            image="HelpIcon.USER_GUIDE">AT User Guide\n');
for m=atchapters()
    fprintf(fid,'            <tocitem target="%s.html">%s</tocitem>\n',m.id,m.title);
    mname = fullfile('m',m.id+".m");
    gid=fopen(mname,'wt');
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
end
