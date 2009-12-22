function atupdate(URLDIR, varargin)
%ATUPDATE syncronizes a local AT installation with SSRL most recent version
%ATUPDATE(URLDIR,OPTION)
% URLDIR    - 'http://www.slac.stanford.edu/~terebilo/public_html/at/code/'
% Options:
% 'info' (default) - Compare local with web and display differences,  no download
%
% 'backup','fullinstall' - not implemented

ATROOTDIR = atroot;

if URLDIR(end)~='/'
    URLDIR = strcat(URLDIR,'/');
end

ATZIPFILE = 'atzip.zip';
ATUPDATEFILE = 'atupdate.txt';



% Make a list of files and updates on SLAC server
S = urlread(strcat(URLDIR,ATUPDATEFILE));

linesep = strfind(S,sprintf('\n'));
timesep = strfind(S,sprintf('\t'));
nlines = length(linesep);

remotefilenames = cell(nlines,1);
remoteupdatetimes = cell(nlines,1);

remotefilenames{1} = S(1:timesep(1)-1);
remoteupdatetimes{1} = S(timesep(1)+1:linesep(1)-1);

for k = 2:length(linesep)
    remotefilenames{k} = S(linesep(k-1)+1:timesep(k)-1);
    remoteupdatetimes{k} = S(timesep(k)+1:linesep(k)-1);
end
    

% Make a list of local files and updates 
fid = fopen('atupdate.txt','r');
if fid>=0
    S = char(fread(fid));
    S = S';
    fclose(fid);

    linesep = strfind(S,sprintf('\n'));
    timesep = strfind(S,sprintf('\t'));
    nlines = length(linesep);

    localfilenames = cell(nlines,1);
    localupdatetimes = cell(nlines,1);

    localfilenames{1} = S(1:timesep(1)-1);
    localupdatetimes{1} = S(timesep(1)+1:linesep(1)-1);

    for k = 2:length(linesep)
        localfilenames{k} = S(linesep(k-1)+1:timesep(k)-1);
        localupdatetimes{k} = S(timesep(k)+1:linesep(k)-1);
    end
else
    localfilenames = {};
    localupdatetimes = {};
end



% New files added on the server since local update
[newfilelist, diffindex] = setdiff(remotefilenames,localfilenames);

if isempty(newfilelist)
    disp('No new files');
end

% Obsolete file - disappeared from the server or moved to a different
% directory or name is now capitlized differently
[obsoletefilelist, diffindex] = setdiff(localfilenames,remotefilenames);

%  Names and directories match
[matchingfilelist, CA, CB] = intersect(localfilenames,remotefilenames);
%  Find feles that more recent on the server
updateindex = find(datenum(remoteupdatetimes(CB))>datenum(localupdatetimes(CA)));
updatefilelist = remotefilenames(CB(updateindex));
up2datefilelist = setdiff(matchingfilelist,updatefilelist);


figure
barh([length(obsoletefilelist), length(newfilelist), length(updatefilelist), length(up2datefilelist)]);
set(gca,'YTickLabel',{'Obsolete Local' , 'New on server', 'Updated on Server', 'Up to date '});
xlabel('Number of files')
title('Comaprison: Your local AT vs. newest AT');


figure
hist(datenum(remoteupdatetimes),30);
datetick('x','mmmyy')
ylabel('Number of  files')
title('File dates of the newest available AT')

return




