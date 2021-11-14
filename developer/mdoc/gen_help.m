function gen_help()
%GEN_HELP	Build the "help" infrastructure
%
%GEN_HELP()     builds the at.m file used by the "help" command

helpfile = fullfile(atroot,'at.m');
atver=ver(atroot);
[fid,fmess]=fopen(helpfile,'wt');
if fid < 0
    error('AT:help','%s: %s',helpfile,fmess);
end

fprintf(fid,'%% %s\n', atver.Name);
fprintf(fid,'%% Version %s %s %s\n%%\n', atver.Version,atver.Release,atver.Date);
fprintf(fid,'%% The Accelerator Toolbox was originally created by Andrei Terebilo.\n');
fprintf(fid,'%% Development is now continued by a multi-laboratory collaboration, %s\n%%\n',...
    weblink('https://github.com/atcollab','atcollab'));
fprintf(fid,'%% %s\n%%\n', weblink('https://atcollab.github.io/at/','AT web site'));
fprintf(fid,'%% The integrators used for tracking particles may be compiled with\n');
fprintf(fid,'%% the command atmexall.\n%%\n');
fprintf(fid,'%% For getting started, one may look at the examples in atmat/atdemos.\n');
fprintf(fid,'%% Example lattice files are located in machine_data.\n');

for m=atchapters()
    fprintf(fid,'%%\n%%*%s*\n',m.title);
    mloop(fid,m.contents);
end

fclose(fid);

    function mloop(fid,mlist)
        bef_fun = true;
        for item=mlist
            if startsWith(item,"-")     % Section header
                fprintf(fid,'%%\n%%   %s\n', eraseBetween(item,1,1));
                bef_fun = false;
            elseif startsWith(item,"0") % Plain text
                fprintf(fid,'%%       %s\n',eraseBetween(item,1,1));
                bef_fun = true;
            else                        % AT function name
                if bef_fun, fprintf(fid,'%%\n'); end
                try
                    h1=h1_line(which(item));
                    line=sprintf('      %-24s - %s',h1.name,h1.h1);
                    fprintf(fid,'%% %s\n',replace(line,h1.name,helplink(h1.name)));
                catch err
                    disp(err.message)
                end
                bef_fun = false;
            end
        end
    end
    function link=helplink(fname)
        link=sprintf('<a href="matlab:help %s">%s</a>', fname,fname);
    end
    function link=weblink(addr,label)
        link=sprintf('<a href="matlab:web(''%s'')">%s</a>',addr,label);
    end
end
