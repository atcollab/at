function gen_contents()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

helpfile = fullfile(atroot,'Contents.m');
[fid,fmess]=fopen(helpfile,'wt');
if fid < 0
    error('AT:help','%s: %s',helpfile,fmess);
end

fprintf(fid,'%% Accelerator toolbox\n');
fprintf(fid,'%% Version 2.0 (atcollab) 01-Aug-2021\n%%\n');
fprintf(fid,'%% The Accelerator Toolbox was originally created by Andrei Terebilo.\n');
fprintf(fid,'%% Development is now continued by a multi-laboratory collaboration, atcollab\n%%\n');
fprintf(fid,'%% The integrators used for tracking particles may be compiled with\n');
fprintf(fid,'%% the command atmexall.\n%%\n');
fprintf(fid,'%% For getting started, one may look at the examples in atmat/atdemos.\n');
fprintf(fid,'%% Example lattice files are located in machine_data.\n');

desc=atchapters();
for m=desc
    fprintf(fid,'%%\n%% *%s*\n%%\n',m.title);
    mloop(fid,m.list);
end
fclose(fid);

    function mloop(fid,mlist)
        for item=mlist
            if startsWith(item,"-")
                fprintf(fid,'%%      %s\n', eraseBetween(item,1,1));
            elseif startsWith(item,"0")
                fprintf(fid,'%% %s\n',eraseBetween(item,1,1));
            else
                try
                    h1=h1_line(which(item));
                    line=sprintf('%-30s - %s',h1.name,h1.h1);
                    % link=sprintf('<matlab:help(''%s'') %s>', h1.name,h1.name);
                    % fprintf(fid,'%% %s\n',replace(line,h1.name,link));
                    fprintf(fid,'%% %s\n',line);
                catch err
                    disp(err.message)
                end
            end
        end
    end
end
