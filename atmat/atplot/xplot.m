function varargout = xplot(ring,is6d,ax,varargin)
%XPLOT  Private function used by atplot and atbaseplot

npts=400; % number of point

% Get options
[synopt,varargs]=getoption(varargin,'synopt',true);
[leftargs,varargs]=getoption(varargs,'leftargs',{});
[rightargs,varargs]=getoption(varargs,'rightargs',{});
[KeepAxis,varargs]=getflag(varargs,'KeepAxis');
[dpargs,varargs]=getoption(varargs,{'dp','df','dct'});
[srange,varargs]=getargs(varargs,[0 inf],'check',@(x) isnumeric(x) && numel(x)==2);

props = atGetRingProperties(ring);
curve.periodicity=props.Periodicity;
elt0=length(ring);

% Select the plotting range
el1=1;
el2=elt0+1;
s0=findspos(ring,el1:el2);
curve.length=s0(end);
smin=max(srange(1), 0);
smax=min(srange(2), curve.length);
els=find(smin>s0,1,'last');
if ~isempty(els), el1=els; end
els=find(s0>smax,1,'first');
if ~isempty(els), el2=els; end

% select the plotting function
[plotfun,varargs]=getargs(varargs, @plotbetadisp,'check', @(x) isa(x,'function_handle'));
[plotargs,varargs]=getargs(varargs,{},'check', @iscell);

dpp=getoption(dpargs,'dp',nan);
if isfinite(dpp) || is6d
    curve.dpp = dpp;
else
    curve.dpp=0;
end

% Split the ring
elmlength=findspos(ring(el1:el2-1),el2-el1+1)/npts;
r2=cellfun(@splitelem,ring(el1:el2-1),'UniformOutput',false);
ring=cat(1,ring(1:el1-1),r2{:},ring(el2:elt0));
plrange=el1:el2+length(ring)-elt0;

[s,outp]=plotfun(ring,curve.dpp,plotargs{:},dpargs{:});
if numel(outp) >= 2
    % plotyy kept instead of yyaxis for octave compatibility...
    [ax2,curve.left,curve.right]=plotyy(ax,...
        s(plrange),outp(1).values(plrange,:),...
        s(plrange),outp(2).values(plrange,:));
    set(ax2(2),'XTick',[],'YColor',get(ax2(1),'YColor'),rightargs{:});
    ylabel(ax2(1),outp(1).axislabel);
    ylabel(ax2(2),outp(2).axislabel);
    linkaxes([ax2(1) ax2(2)],'x');% allows zoom on both right and left plots
elseif numel(outp) == 1
    curve.left=plot(ax,s(plrange),outp(1).values(plrange,:));
    curve.right=[];
    ylabel(ax,outp(1).axislabel);
else
    curve.left=[];
    curve.right=[];
    set(ax,'YLim',[0 1]);
end
set(ax,'XLim',[smin smax],'XGrid','on','YGrid','on',leftargs{:});
xlabel(ax,'s [m]');
if synopt
    curve.lattice=atplotsyn(ax,ring0,varargs{:});  % Plot lattice elements
end
lines=[curve.left;curve.right];
if ~isempty(lines)
    legend(lines,[outp.labels]);
end
if nargout>0, varargout={curve}; end

    function newelems=splitelem(elem)
        if isfield(elem,'Length') && elem.Length > 0
            nslices=ceil(elem.Length/elmlength);
            if ~KeepAxis
                newelems=atdivelem(elem,ones(1,nslices)./nslices);
            else
                newelems=atdivelem(elem,ones(1,nslices)./nslices,'KeepAxis');
            end
        else
            newelems={elem};
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
