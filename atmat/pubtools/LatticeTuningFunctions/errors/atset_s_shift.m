function rerr=atset_s_shift(r,pos,DS)
%
% implements DS longitudinal position drift
% by changing drifts at the sides of the
% elements defined by pos in r
%
% for dipoles the T2(1) field is also changed and the the out DS is
% modified:
% T2(1)=DS*sin(bendignangle)
% DSout=DS*cos(bendignangle)
%
% pos and DS must be the same size
%
%see also: atsetshift atsettilt atsettiltdipole

if length(pos)~=length(DS)
    error('pos and DS must be the same size');
end

rerr=r;

%find dipoles
dipind=findcells(r(pos),'BendingAngle');
posmag=pos;
dsmag=DS;
posdip=pos(dipind);
dsdip=DS(dipind);
posmag(dipind)=[];
dsmag(dipind)=[];

%% STRAIGTH MAGNETS

% find first drift before each element in pos
% and  first drift after  each element in pos
if ~isempty(posmag)
    [driftUP,driftDO]=finddriftsaroundpos(r,posmag);
    
    % shorten drift up
    LDUP0=atgetfieldvalues(r,driftUP,'Length',{1,1});
    rerr=atsetfieldvalues(rerr,driftUP,'Length',LDUP0-dsmag');
    % lengthen drift down
    LDDO0=atgetfieldvalues(r,driftDO,'Length',{1,1});
    rerr=atsetfieldvalues(rerr,driftDO,'Length',LDDO0+dsmag');
    rerr=atsetfieldvalues(rerr,posmag,'DeltaS',dsmag');
    
end
if ~isempty(posdip)
    %% DIPOLES
    
    % find first drift before each element in pos
    % and  first drift after  each element in pos
    [driftUP,driftDO]=finddriftsaroundpos(r,posdip);
    theta=atgetfieldvalues(rerr,posdip,'BendingAngle',{1,1});
    
    %if dipoles have the same MagNum, move each part and sum effect on T2
    %coordinate change. DS assigned to first slice is assumed as DS of the
    %whole magnet.
    
    maggr=getMagGroupsFromMagNum(r(posdip));
    
    for imaggr=1:length(maggr)
        dipind=maggr{imaggr}; % sliced dipole indexes in posdip
        
        % shorten drift up for first dipole in group
        LDUP0=atgetfieldvalues(rerr,driftUP(dipind(1)),'Length',{1,1});
        rerr=atsetfieldvalues(rerr,driftUP(dipind(1)),'Length',LDUP0-dsdip(dipind(1))');
        
        dsout=dsdip((dipind(1)))';
        dt2out=dsdip((dipind(1)))';
        
        for iiind=1:length(dipind)
            dsout=dsout.*cos(theta(dipind(iiind)));
            dt2out=dt2out.*sin(theta(dipind(iiind)));
        end
        
        try
            dt2out0=atgetfieldvalues(rerr,posdip(dipind(end)),'T2',{1,1}); % add to existing if exists
        catch
            dt2out0=zeros(dt2out);
        end
        
        % lengthen drift down FOR last DIPOLE in group. ALSO T2 changes!
        LDDO0=atgetfieldvalues(rerr,driftDO(dipind(end)),'Length',{1,1});
        rerr=atsetfieldvalues(rerr,driftDO(dipind(end)),'Length',LDDO0+dsout);%+dsdip(dipind(1)));%
        rerr=atsetfieldvalues(rerr,posdip(dipind(end)),'T2',{1,1},dt2out0-dt2out); %
        
        rerr=atsetfieldvalues(rerr,posdip(dipind),'DeltaS',dsdip(dipind(1))');
        rerr=atsetfieldvalues(rerr,posdip(dipind),'DeltaST2',...
            [zeros(1,length(dipind)-1) dt2out]);
        
    end
        
end

return


function [dup,ddo]=finddriftsaroundpos(r,pos)
dup=nan(size(pos));
ddo=nan(size(pos));

for indpos=1:length(pos)
    
    i=pos(indpos);
    
    NEL=length(r);
    while ~strcmp(r{i}.Class,'Drift')
        if i<NEL
            i=i+1;
        else
            i=1;
        end
    end
    dup(indpos)=i;
    
    i=pos(indpos);
    while ~strcmp(r{i}.Class,'Drift')
        if i>1
            i=i-1;
        else
            i=NEL;
        end
    end
    ddo(indpos)=i;
end

return


function a=getmagnumdipole(r,ind)

try
    a=r{ind}.MagNum;
catch
    a=NaN;
end

return