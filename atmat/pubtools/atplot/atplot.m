function curve = atplot(varargin)
%ATPLOT Plots optical functions
%
%ATPLOT                 Plots THERING in the current axes
%
%ATPLOT(RING)           Plots the lattice specified by RING
%
%ATPLOT(AX,RING)        Plots in the axes specified by AX
%
%ATPLOT(AX,RING,DPP)    Plots at momentum deviation DPP
%
%ATPLOT(...,[SMIN SMAX])  Zoom on the specified range
%
%ATPLOT(...,@PLOTFUNCTION,args...)
%	Allows for a user supplied function providing the values to be plotted
%   PLOTFUNCTION must be of form:
%   PLOTDATA=PLOTFUNCTION(LINDATA,RING,DPP,args...)
%
%       PLOTDATA is a structure array,
%       PLOTDATA(1) describes data for the left axis
%           PLOTDATA(1).VALUES: data to be plotted length(ring)+1 X nbcurve
%           PLOTDATA(1).LABELS: curve labela, cell array, 1 X nbcurve
%           PLOTDATA(1).AXISLABEL: string
%       PLOTDATA(2) optional, describes data for the left axis
%
%   The default function displayed below as an example plots beta functions
%	and dispersion
%
% function plotdata=defaultplot(lindata,ring,dpp,varargin)
% beta=cat(1,lindata.beta);                     % left axis
% plotdata(1).values=beta;
% plotdata(1).labels={'\beta_x','\beta_z'};
% plotdata(1).axislabel='\beta [m]';
% dispersion=cat(2,lindata.Dispersion)';        % right axis
% plotdata(2).values=dispersion(:,1);
% plotdata(2).labels={'\eta_x'};
% plotdata(2).axislabel='dispersion [m]';
% end
%
%CURVE=ATPLOT(...) Returns handles to some objects:
%   CURVE.LEFT      Handles to the left axis plots
%   CURVE.RIGHT     Handles to the right axis plots
%   CURVE.LATTICE	Handles to the Element patches
%   CURVE.COMMENT	Handles to the Comment text
%
%See also: atplot>defaultplot
global THERING

narg=1;
plotfun=@defaultplot;
% Select axes for the plot
if narg<=length(varargin) && isscalar(varargin{narg}) && ishandle(varargin{narg});
    ax=varargin{narg};
    narg=narg+1;
else
    ax=gca;
end
% Select the lattice
if narg<=length(varargin) && iscell(varargin{narg});
    [elt0,nperiods,ring0]=get1cell(varargin{narg});
    narg=narg+1;
else
    [elt0,nperiods,ring0]=get1cell(THERING);
end

% Select the momentum deviation
if narg<=length(varargin) && isscalar(varargin{narg}) && isnumeric(varargin{narg});
    dpp=varargin{narg};
    narg=narg+1;
else
    dpp=0;
end
s0=findspos(ring0,1:elt0+1);
srange=s0([1 end]);     %default plot range
el1=1;
el2=elt0+1;

firstarg=nargin+1;
for iarg=narg:nargin
    % explicit plot range
    if isnumeric(varargin{iarg}) && (numel(varargin{iarg})==2)
        srange=varargin{iarg};
        els=find(srange(1)>s0,1,'last');
        if ~isempty(els), el1=els; end
        els=find(s0>srange(2),1,'first');
        if ~isempty(els), el2=els; end
    elseif isa(varargin{iarg},'function_handle')
        plotfun=varargin{iarg};
        firstarg=iarg+1;
        break
    else
        firstarg=iarg;
        break
    end
end

%ring=[ring0((1:el1-1)');atslice(ring0(el1:el2-1),250);ring0((el2:elt0)')];
ring=[ring0(1:el1-1,1);atslice(ring0(el1:el2-1,1),250);ring0(el2:elt0,1)];
elt=length(ring);
plrange=el1:el2+elt-elt0;
[lindata,tuneper,chrom]=atlinopt(ring,dpp,1:elt+1); %#ok<NASGU,ASGLU>
s=cat(1,lindata.SPos);

%axclean(ax);
set(ax,'Position',[.13 .11 .775 .775],'FontSize',12);
outp=plotfun(lindata,ring,dpp,varargin{firstarg:nargin});
labels=[outp.labels];
if numel(outp) >= 2
    [ax2,h1,h2]=plotyy(ax,s(plrange),outp(1).values(plrange,:),s(plrange),outp(2).values(plrange,:));
    set(ax2(2),'XTick',[],'FontSize',12);
    ylabel(ax2(1),outp(1).axislabel);
    ylabel(ax2(2),outp(2).axislabel);
else
    h1=plot(ax,s(plrange),outp(1).values(plrange,:));
    h2=[];
    ax2=ax;
    ylabel(ax,outp(1).axislabel);
end
curve.left=h1;
curve.right=h2;
curve.lattice=atplotsyn(ax2(1),ring0);  % Plot lattice elements
lts=get(ax2(1),'Children');             % Put them in the background
set(ax2(1),'Children',[lts(4:end);lts(1:3)]);
set([h1;h2],'LineWidth',1);
set(ax2,'XLim',srange);
xlabel(ax2(1),'s [m]');                 % Add labels
legend([h1;h2],labels{:});
grid on
% axes(ax);
tuneper=lindata(end).mu/2/pi;
tunes=nperiods*tuneper;
circ=nperiods*s0(end);
if nperiods > 1, plural='s'; else plural=''; end
line1=sprintf('\\nu_x=%8.3f      \\deltap/p=%.3f%i %s',tunes(1),dpp);
line2=sprintf('\\nu_z=%8.3f      %2i %s, C=%10.3f',tunes(2),nperiods,...
    ['period' plural],circ);
curve.comment=text(-0.14,1.12,{line1;line2},'Units','normalized',...
    'VerticalAlignment','top','FontSize',12);

    function axclean(ax)
        cla(ax,'reset');
        appd=fieldnames(getappdata(ax));
        for as=appd'
            rmappdata(ax,as{1});
        end
        fig=get(ax,'Parent');
        axs=get(fig,'Children');
        for testax=axs'
            testax
            getappdata(testax)
        end
    end

    function [cellsize,np,cell]=get1cell(ring)
        [cellsize,np]=size(ring);
        cell=ring(:,1);
        params=atgetcells(cell,'Class','RingParam');
        if any(params)
            np=ring{find(params,1)}.Periodicity;
        end
    end
end

function plotdata=defaultplot(lindata,ring,dpp,varargin) %#ok<INUSD>
%DEFAULTPLOT    Default plotting function for ATPLOT
%
%Plots beta-functions on left axis and dispersion on right axis

beta=cat(1,lindata.beta);
plotdata(1).values=beta;
plotdata(1).labels={'\beta_x','\beta_z'};
plotdata(1).axislabel='\beta [m]';
dispersion=cat(2,lindata.Dispersion)';
plotdata(2).values=dispersion(:,1);
plotdata(2).labels={'\eta_x'};
plotdata(2).axislabel='dispersion [m]';
end
