function [superp,periods]=atreadbeta(filename,cavipass,bendpass,quadpass)
%ATREADBETA				reads a BETA file
%
%ring=ATREADBETA(fname,cavipass,bendpass,quadpass,multipass)
%
%FILENAME:	BETA file
%CAVIPASS:	pass method for cavities (default IdentityPass)
%BENDPASS:	pass method for dipoles (default BndMPoleSymplectic4E2Pass)
%QUADPASS:	pass method for quadrupoles (default QuadMPoleFringePass)
%MULTIPASS:	pass method for sextupoles (default StrMPoleSymplectic4Pass)
%
%[superp,periods]=ATREADBETA(fname,cavipass,bendpass,quadpass)
%		returns only one superperiod and the number of superperiods
%
%See also: ATX, ATLINOPT, ATMODUL

global GLOBVAL
persistent fpath

if isempty(fpath), fpath=getenv('DBETA'); end

if nargin < 5, multipass='StrMPoleSymplectic4Pass'; end
if nargin < 4, quadpass='QuadMPoleFringePass'; end
if nargin < 3, bendpass='BndMPoleSymplectic4E2Pass'; end
if nargin < 2, cavipass='IdentityPass'; end
if nargin < 1, filename=''; end

if isempty(filename)
    [fname,fpath]=uigetfile('*.str','BETA structure',[fpath filesep]);
    if ~ischar(fname), error('ReadBeta:NoFile','No file selected'); end
    filename=fullfile(fpath,fname);
end

fid=fopen(filename,'rt');
if fid < 0, error('ReadBeta:NoFile','Cannot open file %s',filename); end

betadelim(fid,'LIST OF ELEMENTS');
GLOBVAL.E0=1E9;
line=fgetl(fid);
nb_elems=sscanf(line,'%d');
elemtable=struct();
cavilist={};
for el=1:nb_elems
    nextelem=readelem(fid,cavipass,bendpass,quadpass,multipass);
    try
        elemtable.(nextelem.FamName)=nextelem;
    catch %#ok<CTCH>
        nextelem.FamName=['x' nextelem.FamName];
        elemtable.(nextelem.FamName)=nextelem;
    end
    if strcmp(nextelem.BetaCode,'CA')
        cavilist{end+1}=nextelem; %#ok<AGROW>
    end
end
if isempty(cavilist)% Prepare a default cavity
    cavilist{1}=rfcavity('RFCAV',0,0,0,1,cavipass);
end
eledict=fieldnames(elemtable);
disp(['Elements processed (' num2str(nb_elems) ' elements)']);

betadelim(fid,'STRUCTURE');
line=fgetl(fid);
nb_stru=sscanf(line,'%d');
superp=cell(nb_stru,1);

dipelem=[];         % Pending dipole element (waiting for exit face)
anglein=0;          % Current dipole face angles
angleout=0;
lff=0;              % Current dipole fringe field extension
displ=zeros(1,3);   % Current misalignment vector
srot=0;             % Current element rotation

id_stru=0;
for el=1:nb_stru
    elcode=fscanf(fid,'%s',1);              % Select element in the table
    try
        elnum=str2double(elcode);
        if isfinite(elnum)
            elcode=eledict{elnum};
        end
        nextelem=elemtable.(elcode);
    catch %#ok<CTCH>
        error('ReadBeta:BadElem',['Cannot identify element ' elcode]);
    end
    switch nextelem.BetaCode                % Process the element
        case 'CO'
            if isempty(dipelem)             % Entrance face
                anglein=nextelem.Angle;
                lff=nextelem.Lff;
            else
                angleout=nextelem.Angle;    % Exit face
                id_stru=id_stru+1;          % create immediately in case of 2 adjacent CO elems
                superp{id_stru}=atelem(dipelem,'EntranceAngle',anglein,...
                    'ExitAngle',angleout,'FullGap',0,'FringeInt',lff);
                anglein=0;
                angleout=0;
                lff=0;
                dipelem=[];
            end
        case 'RO'
            srot=srot+nextelem.Srot;
        case 'DE'
            displ=displ+nextelem.Displacement;
        otherwise
            if ~isempty(dipelem)
                id_stru=id_stru+1;
                superp{id_stru}=atelem(dipelem,'EntranceAngle',anglein,...
                    'ExitAngle',angleout,'FullGap',0,'FringeInt',lff);
                anglein=0;
                angleout=0;
                lff=0;
                dipelem=[];
            end
            if srot ~= 0
                if isfield(nextelem,'R1')
                    srollmat=mkSRotationMatrix(srot);
                    nextelem.R1=srollmat;
                    nextelem.R2=srollmat';
                end
            end
            if max(abs(displ)) > 0
                if isfield(nextelem,'T1')
                    nextelem.T1([1 3 5])=displ;
                    nextelem.T2([1 3 5])=-displ;
                end
            end
            if strcmp(nextelem.BetaCode,'DI')
                dipelem=nextelem;
            else
                id_stru=id_stru+1;
                superp{id_stru}=nextelem;
            end
    end
end
superp(id_stru+1:end)=[];
disp(['Structure processed (' num2str(nb_stru) ' elements)']);

nper=fscanf(fid,'%d',1);
fclose(fid);

cavities=findcells(superp,'BetaCode','CA');
if isempty(cavities)		% add implicit cavity if necessary
    superp{end+1}=cavilist{1};
    cavities=length(superp);
end
superp=setcellstruct(superp,'Energy',true(size(superp)),GLOBVAL.E0);

if nargout >= 2
    periods=nper;
    for i=cavities			% set cavity frequency
        superp{i}=tunecavity(superp{i},9E6/16/length(cavities),...
            findspos(superp,id_stru+1),1,nper);
    end
else
    for i=cavities			% set cavity frequency
        superp{i}=tunecavity(superp{i},9E6/16/length(cavities),...
            findspos(superp,id_stru+1),nper,nper);
    end
    if nper > 1
        superp=repmat(superp,1,nper);
    end
end


evalin('base','global GLOBVAL');

function cav=tunecavity(cav,V,clength,ncell,ntot)
frev=2.997924E8/clength;
if cav.Voltage == 0, cav.Voltage=V; end
if cav.HarmNumber > 1
    harm=ceil(cav.HarmNumber/ntot);
else
    harm=round(1.66678*clength);
end
cav.Frequency=frev*harm;
cav.HarmNumber=ncell*harm;

function newelem=readelem(fid,cavipass,bendpass,quadpass,multipass)

global GLOBVAL

line=fgetl(fid);
next=1;
[elname,count,errmess,nl]=sscanf(line(next:end),'%s',1); %#ok<ASGLU>
next=next+nl;
[code,count,errmess,nl]=sscanf(line(next:end),'%s',1); %#ok<ASGLU>
next=next+nl;
params=sscanf(line(next:end),'%f')';
switch (code)
    case 'SD'
        newelem=atdrift(elname,params(1));
    case 'QP'
        newelem=atquadrupole(elname,params(1),params(2),quadpass);
    case 'DE'
        newelem=atelem(atmarker(elname),'Displacement',params(1:3));
    case 'RO'
        newelem=atelem(atmarker(elname),'Srot',params(1));
    case 'CO'
        newelem=atelem(atmarker(elname),'Angle',params(1),'Lff',params(3));
    case 'DI'
        strength=-params(3)/params(2)/params(2);
        newelem=atsbend(elname,params(1)*params(2),params(1),strength,bendpass);
    case 'SX'
        if params(1) < 0.001
            code='LD3';
            newelem=atthinmultipole(elname,[],[0 0 params(1)*params(2)]);
        else
            newelem=atsextupole(elname,params(1),params(2),multipass);
        end
    case 'LD'
        order=params(2)/2;
        polb=[];
        polb(1,order)=params(1);
        code=[code int2str(order)];
        newelem=atthinmultipole(elname,[],polb);
    case 'LT'
        order=params(2)/2;
        pola=[];
        pola(1,order)=params(1);
        code=[code int2str(order)];
        newelem=atthinmultipole(elname,pola,[]);
    case 'KI'
        if params(3) > 0, code='CHV'; end
        newelem=atelem('corrector','FamName',elname,'KickAngle',[params(1) params(2)]);
        newelem.PassMethod='IdentityPass';
    case 'CA'
        GLOBVAL.E0=params(3);
        newelem=atelem(atmarker(elname,cavipass),'Length',0,...
            'Voltage',abs(params(1)),'Frequency',0,'HarmNumber',params(2));
    otherwise
        newelem=atmarker(elname);
end
newelem.BetaCode=code;

function betadelim(fid,code)
while true
    while true
        line=fgetl(fid);
        if ~ischar(line), error('ReadBeta:EndOfFile','Encountered unexpected end of file.'); end
        if ~isempty(strfind(line,'***')), break; end
    end
    if ~isempty(strfind(line,code)), break; end
end
