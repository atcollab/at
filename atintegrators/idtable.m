function [z L] = idtable(fname, Nslice, filename, Energy, method)
% idtable(fname, Nslice, filename, Energy, method)
%
% FamName	family name
% Nslice	number of slices (1 means the wiggler is represented by a
%           single kick in the center of the device).
% filename	name of file with wiggler tracking tables.
% Energy    Energy of the machine, needed for scaling
% method    name of the function to use for tracking. Use 'WigTablePass'
%
% The tracking table is described in
% P. Elleaume, "A new approach to the electron beam dynamics in undulators
% and wigglers", EPAC92.
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% 13-09-2007:  Created by M. Munoz, based in J. Safranek code.
% 17-11-2008:  Modificated by Z.Martí
% 20-07-2017:  I. Martin
%---------------------------------------------------------------------------

brho = 1e9/299792458*Energy;

fid = fopen(filename, 'r');
data = textscan(fid,'%s','delimiter','\n');
data = data{1};
fclose(fid);

% parse the input data for key parameters
startloc = [];
for i = 1:length(data);
    if strfind(data{i},'Undulator'); L = str2num(data{i+1}); end
    if strfind(data{i},'Horizontal'); if strfind(data{i},'Points'); nx = str2num(data{i+1});  end; end
    if strfind(data{i},'Vertical'  ); if strfind(data{i},'Points'); ny = str2num(data{i+1});  end; end
    if strfind(data{i},'START'); startloc(end+1) = i+1; end
end

% initialise arrays
xkick  = zeros(ny,nx);
ykick  = zeros(ny,nx);
skick  = zeros(ny,nx);
xkick1 = zeros(ny,nx);
ykick1 = zeros(ny,nx);

block = 1;
xtable = str2num(data{startloc(block)});
for i = startloc(block) + 1:startloc(block) + ny;
    temp = str2num(data{i});
    ytable(i-startloc(block)) = temp(1);
    xkick(i-startloc(block),1:nx) = temp(2:end)/brho/brho;
end
block = 2;
for i = startloc(block) + 1:startloc(block) + ny;
    temp = str2num(data{i});
    ykick(i-startloc(block),1:nx) = temp(2:end)/brho/brho;
end

% add in the optional fields
if length(startloc)==3
    % if there are 3 blocks of data, the third is skick
    block = 3;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        skick(i-startloc(block),1:nx) = temp(2:end)/brho/brho;
    end
elseif length(startloc)==4
    % if there are 4 blocks of data, the third is xkick1 and the fourth ykick1
    block = 3;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        xkick1(i-startloc(block),1:nx) = -1*temp(2:end)/brho;
    end
    block = 4;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        ykick1(i-startloc(block),1:nx) = -1*temp(2:end)/brho;
    end
elseif length(startloc)==5
    % if there are 5 blocks of data, the third is skick, the fourth is xkick1 and the fifth is ykick1
    block = 3;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        skick(i-startloc(block),1:nx) = temp(2:end)/brho/brho;
    end
    block = 4;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        xkick1(i-startloc(block),1:nx) = -1*temp(2:end)/brho;
    end
    block = 5;
    for i = startloc(block) + 1:startloc(block) + ny;
        temp = str2num(data{i});
        ykick1(i-startloc(block),1:nx) = -1*temp(2:end)/brho;
    end
end

% Sort arrays in ascending order (needed for "IdTablePass.c")
[ytable indy]=sort(ytable);
[xtable indx]=sort(xtable);
xtable = xtable(:);
ytable = ytable(:);
xkick=xkick(indy,indx);
ykick=ykick(indy,indx);


ElemData.FamName = fname;  % add check for identical family names
ElemData.Length = L;
ElemData.Nslice = Nslice;
ElemData.xtable = xtable;
ElemData.ytable = ytable;
ElemData.xkick = xkick;
ElemData.ykick = ykick;
ElemData.skick = skick;
ElemData.xkick1 = xkick1;
ElemData.ykick1 = ykick1;
ElemData.PassMethod = method;
ElemData.R1 = diag(ones(6,1));
ElemData.R2 = diag(ones(6,1));
ElemData.T1 = zeros(1,6);
ElemData.T2 = zeros(1,6);



global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;
