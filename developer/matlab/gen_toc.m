function gen_toc()
%GEN_TOC	Build the HTML files used by the Matlab help browser

[devdir,~, ~]=fileparts(mfilename('fullpath'));
docdir = fullfile(atroot,'..','docs','matlab');
tocfile = fullfile(docdir,'helptoc.xml');
[fid,fmess,fmid] = copyfile(fullfile(devdir,'helptoc.xml'),tocfile);
if fid <= 0
    error(fmid,'%s: %s',tocfile,fmess);
end
fid=openmfile(tocfile,'at');

% User guide
ugname=fullfile('m','ugsummary.m');
hid=openmfile(ugname,'wt');
fprintf(fid,'        <tocitem target="ugsummary.html"\n');
fprintf(fid,'            image="HelpIcon.USER_GUIDE">AT User Guide\n');
fprintf(hid,'%%%% AT User Guide\n%%\n%%%%\n');
%   Loop on UG chapters
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
fclose(hid);

% How to...
mlxloop('howtos','How to…',@howtochapters);

% Release notes
mlxloop('release_notes','Release notes');

% Web site
fprintf(fid,'        <tocitem target="https://atcollab.github.io/at/" \n');
fprintf(fid,'                 image="$toolbox/matlab/icons/webicon.gif">\n');
fprintf(fid,'        AT Web Site\n');
fprintf(fid,'        </tocitem>\n');

fprintf(fid,'    </tocitem>\n');
fprintf(fid,'</toc>\n');
fclose(fid);

publish(ugname,'evalCode',false,'outputDir',docdir);

% Publish custom files
for dd=reshape(dir(fullfile(devdir,'mlx')),1,[])
    [~,nn,xx]=fileparts(dd.name);
    if strcmp(xx,'.mlx')
        export(fullfile(dd.folder,dd.name),fullfile(docdir,strcat(nn,'.html')));
    end
end

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

    function mlxloop(secdir,secname,chapfun)
        if nargin<3, chapfun=@lst; end
        dirname=fullfile(devdir,secdir);
        sumname=fullfile(devdir,'m',[secdir '.m']);
        sumid=openmfile(sumname,'wt');
        fprintf(fid,'        <tocitem target="%s" \n', [secdir '.html']);
        fprintf(fid,'                 image="$toolbox/matlab/icons/webicon.gif">\n');
        fprintf(fid,'        %s\n',secname);
        fprintf(sumid,'%%%% %s\n%% \n%%%%\n', secname);
        for mm=chapfun()
            target=fullfile(secdir,mm.id+".html");
            export(fullfile(secdir,mm.id+".mlx"),fullfile(docdir,target));
            fprintf(sumid,'%% <matlab:web(fullfile(docroot,''3ptoolbox'',''atacceleratortoolbox'',''doc'',''%s'')) %s>\n%%\n',target,mm.title);
            fprintf(fid,'            <tocitem target="%s">%s</tocitem>\n',target,mm.title);
        end
        fclose(sumid);
        publish(sumname,'evalCode',false,'outputDir',docdir);
        fprintf(fid,'        </tocitem>\n');
        delete(sumname);
        function res=lst()
            vals=reshape(dir(fullfile(dirname,'*.mlx')),1,[]);
            [~,nms,~]=arrayfun(@fileparts,{vals.name},'UniformOutput',false);
            nms=string(sort(nms));
            res=struct('id',nms,'title',nms);
        end
    end

    function fid=openmfile(fpath,mode)
        [fid,mess]=fopen(fpath,mode);
        if fid < 0
            error('AT:doc','Cannot open %s: %s',fpath,mess);
        end
    end

end
