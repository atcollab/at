function setversion(version)
%SETVERSION Set the version of AT
%
%SERVERSION(VERSION)

versfile=fullfile(atroot,'Contents.m');
tmpfile=fullfile(tempdir,'Contents.m');
fout=fopen(tmpfile,'wt');
fin=fopen(versfile,'rt');
line=fgetl(fin); %#ok<NASGU> 
line=fgetl(fin); %#ok<NASGU> 
line=fgetl(fin);
fprintf(fout,'%% Accelerator Toolbox\n');
fprintf(fout,'%% Version %s (atcollab) %s\n',version,datestr(datetime,1));
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