function curve = atplot(varargin)
%ATPLOT Plots optical functions
%
%ATPLOT             Plots THERING in the current axes
%
%ATPLOT(RING)       Plots the lattice specified by RING
%
%ATPLOT(AX,RING)    Plots in the axes specified by AX
%
%ATPLOT(...,[SMIN SMAX])  Zoom on the specified range
%
%CURVE=ATPLOT(...) Returns handles to some objects:
%   CURVE.BETAX    Hor. beta function line
%   CURVE.BETAZ    Vert. beta function line
%   CURVE.ETAX     Hor. dispersion line
%   CURVE.LATTICE  Element patches
%   CURVE.COMMENT  Comment text

global THERING

narg=1;
                        % Select axes for the plot
if narg<=length(varargin) && isscalar(varargin{narg}) && ishandle(varargin{narg});
    ax=varargin{narg};
    narg=narg+1;
else
    ax=gca;
end
                        % Select the lattice
if narg<=length(varargin) && iscell(varargin{narg});
    ring0=varargin{narg};
    narg=narg+1;
else
    ring0=THERING;
end
[elt0,nperiods]=size(ring0);
s0=findspos(ring0,1:elt0+1);
srange=s0([1 end]);     %default plot range
el1=1;
el2=elt0+1;
for iarg=narg:nargin
                        % explicit plot range
    if isnumeric(varargin{iarg}) && (numel(varargin{iarg})==2)
        srange=varargin{iarg};
        el1=find(srange(1)>s0,1,'last');
        el2=find(s0>srange(2),1,'first');
    end
end
%ring=[ring0((1:el1-1)');atslice(ring0(el1:el2-1),250);ring0((el2:elt0)')];
ring=[ring0(1:el1-1,1);atslice(ring0(el1:el2-1),250);ring0(el2:elt0,1)];
elt=length(ring);
plrange=el1:el2+elt-elt0;
[lindata,tunes,chrom]=atlinopt(ring,0,1:elt+1); %#ok<NASGU,ASGLU>
beta=cat(1,lindata.beta);
disp=cat(2,lindata.Dispersion)';
s=cat(1,lindata.SPos);
set(ax,'Position',[.13 .11 .775 .775]);
[ax2,h1,h2]=plotyy(ax,s(plrange),beta(plrange,:),s(plrange),disp(plrange,1));
set([h1;h2],'LineWidth',1);
set(ax2(2),'XTick',[]);
curve.betax=h1(1);
curve.betaz=h1(2);
curve.etax=h2(1);
set(ax2,'XLim',srange);
curve.lattice=atplotsyn(ax2(1),ring0);  % Plot lattice elements
lts=get(ax2(1),'Children');             % Put them in the background
set(ax2(1),'Children',[lts(4:end);lts(1:3)]);
xlabel(ax2(1),'s [m]');                 % Add labels
ylabel(ax2(1),'\beta [m]');
ylabel(ax2(2),'\eta [m]');
legend([h1;h2],'\beta_x','\beta_z','\eta_x');
grid on
% axes(ax);
tuneper=lindata(end).mu/2/pi;
tunes=nperiods*tuneper;
if nperiods > 1, plural='s'; else plural=''; end
line1=sprintf('\\nu_x=%8.3f     %i %s',tunes(1),nperiods,['period' plural]);
line2=sprintf('\\nu_z=%8.3f     C=%10.3f',tunes(2),nperiods*s0(end));
curve.comment=text(-0.14,1.12,{line1;line2},'Units','normalized',...
    'VerticalAlignment','top');
end
