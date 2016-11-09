function [bestinputcoord]=Scan2x2DinCOD(ropen,inCOD,ngridstep,gridsize,figlabel)
%[bestinputcoord]=ScanPosAngle(...
%     ropen,... lattice 
%     inCOD,... initial trajectory [x,xp,y,yp]
%     ngridstep,... grid steps in 1D (51 by default)
%     gridsize,... grid size (5mm by default)
%     figlabel)     figure label name. if not given, no figure
%
% scans initial coordinates for maximum number of turns
%
%see also: ringpass

if nargin<5
figlabel=[];
end
if nargin<4
    gridsize=5e-3;
end
if nargin<3
    ngridstep=51;
end

% maxexcurs=20e-3;
% Xapert=maxexcurs*ones(size(ropen));
% Yapert=maxexcurs*ones(size(ropen));
% ropen=SetPhysicalAperture(ropen,Xapert/2,Yapert/2);

bestinputcoord=inCOD*0;
TotTurn=10;

DX=linspace(-1,1,ngridstep)*gridsize;
DXP=DX/5;
DY=DX;
DYP=DX;

% initial coordinate matrix

[x,xp]=meshgrid(DX,DXP);
[y,yp]=meshgrid(DY,DYP);

incoory=repmat(inCOD,1,length(x(:)))+...
    [repmat(bestinputcoord([1]),length(x(:)),1),...
    repmat(bestinputcoord([2]),length(x(:)),1),...
    y(:),yp(:),zeros(size(xp(:),1),2)]';

[~,lost,NT,lossinfo]=ringpass(ropen,incoory,TotTurn);

lossinfo.element(isnan(lossinfo.element))=0;% particle not lost!
NT(lost==0)=TotTurn;
maxturnsy=(max(NT)+1);
elpasy=lossinfo.element+length(ropen).*NT;

[~,indt]=max(elpasy);
% if more then one point, take closest to center.
guessbest=sort([indt find(lost==0)]);

dist=sqrt(incoory([3],guessbest).^2+incoory([4],guessbest).^2);
[~,indb]=min(dist);

bestinputcoord([3,4])=incoory([3,4],guessbest(indb));


incoorx=repmat(inCOD,1,length(x(:)))+...
    [x(:),xp(:),...
    repmat(bestinputcoord([3]),length(x(:)),1),...
    repmat(bestinputcoord([4]),length(x(:)),1),...
    zeros(size(xp(:),1),2)]';

[~,lost,NT,lossinfo]=ringpass(ropen,incoorx,TotTurn);
lossinfo.element(isnan(lossinfo.element))=0;% particle not lost!
NT(lost==0)=TotTurn;
maxturnsx=(max(NT)+1);

elpasx=lossinfo.element+length(ropen).*NT;

[~,indt]=max(elpasx);
guessbest=sort([indt find(lost==0)]);

dist=sqrt(incoorx([1],guessbest).^2+incoorx([2],guessbest).^2);
[~,indb]=min(dist);

bestinputcoord([1,2])=incoorx([1,2],guessbest(indb));


%%
if ~isempty(figlabel)
    maxturns=max([maxturnsx,maxturnsy]);
    
    figure('Units','normalized','Position',[0.1 0.2 0.6 0.5]);
    subplot(1,2,1);
    surf(x,xp,reshape(elpasx,length(DX),[]));
      hcb=colorbar;
    caxis([0,length(ropen)*maxturns]);
    set(hcb,'YTick',[1:maxturns]*length(ropen));
    set(hcb,'YTickLabel',arrayfun(@(a)num2str(a,'turn %d'),[1:maxturns],'un',0));
    xlabel('x');ylabel('xp');view(2);shading flat;
    %title({['Number of turns vs injection position'],['Initial: ' num2str(inCOD)],['Best: ' num2str(bestinputcoord) ]});
    title({[],[],[],[],[]})
    subplot(1,2,2);
    surf(y,yp,reshape(elpasy,length(DY),[]));
    xlabel('y');ylabel('yp');view(2);shading flat;
     hcb=colorbar;
    caxis([0,length(ropen)*maxturns]);
    set(hcb,'YTick',[1:maxturns]*length(ropen));
    set(hcb,'YTickLabel',arrayfun(@(a)num2str(a,'turn %d'),[1:maxturns],'un',0));
    title({[],[],[],[],[]})
    
    set(gcf,'NextPlot','add');
    axes;
    h = title({['Number of turns vs injection position'],['Initial: ' num2str(inCOD)],['Best: ' num2str(bestinputcoord) ]});
    set(gca,'Visible','off');
    set(h,'Visible','on');

    
    saveas(gca,[figlabel '.fig']);
    try
        export_fig([figlabel '.jpg']);
    catch
    disp('missing export_fig, no jpg')
    end
    
end

return