function pp = atplotsyn(ax,ring,varargin)
%ATPLOTSYN Helper function for ATPLOT
%
%PATCHES=ATPLOTSYN(AX,RING) Plots the magnetic elements found in RING

options=struct(varargin{:});
labs=false(size(ring(:)));
if isfield(options,'labels')
    if islogical(options.labels)
        labs=options.labels(:);
    elseif isnumeric(options.labels)
        labs(options.labels)=true;
    end
end
xlim=get(ax,'XLim');
slim=diff(xlim);
tlim=diff(get(ax,'YLim'));
nlines=length(get(ax,'Children'));
%axes(ax);

sl=findspos(ring(:,1),1:size(ring,1)+1);
ll=diff(sl);
sok=(sl(2:end)>xlim(1)) & (sl(1:end-1)<xlim(2));
sl=sl(sok);
ll=ll(sok);
rok=ring(sok);
labs=labs(sok);
[dipoles,qpoles,spoles,mpoles,bpms]=cellfun(@eltype,rok);

ampl=0.05*tlim;     % Dipoles
[xplot,yplot]=setxpl(sl(dipoles),ll(dipoles),[0;0;1;1;0],ampl*[0;1;1;0;0]);
%[xplot,yplot]=setxpl(sl(dipoles),ll(dipoles),[0;0;1;1;0],ampl*[-1;1;1;-1;-1]);
p1=patch(xplot,yplot,[0.5 0.5 1],'DisplayName','Dipoles');

ampl=0.05*tlim;     % Quadrupoles
foc=reshape(0.4*ampl*sign(atgetfieldvalues(rok(qpoles),'PolynomB',{2})),1,[]);
[xplot,yplot]=setxpl(sl(qpoles),ll(qpoles),[0;0;0.5;1;1;0],ampl*[0;1;1;1;0;0]);
yplot(3,:)=yplot(3,:)+foc;
%[xplot,yplot]=setxpl(sl(qpoles),ll(qpoles),[0;0;0.5;1;1;0.5;0],ampl*[-1;1;1;1;-1;-1;-1]);
%yplot(6,:)=yplot(6,:)-foc;
p2=patch(xplot,yplot,[1 0.5 0.5],'DisplayName','Quadrupoles');

ampl=0.04*tlim;     % Sextupoles
foc=reshape(0.25*ampl*sign(atgetfieldvalues(rok(spoles),'PolynomB',{3})),1,[]);
[xplot,yplot]=setxpl(sl(spoles),ll(spoles),[0;0;0.33;0.66;1;1;0],ampl*[0;1;1;1;1;0;0]);
yplot(3:4,:)=yplot(3:4,:)+foc([1;1],:);
%[xplot,yplot]=setxpl(sl(spoles),ll(spoles),[0;0;0.33;0.66;1;1;0.66;0.33;0],ampl*[-1;1;1;1;1;-1;-1;-1;-1]);
%yplot(7:8,:)=yplot(7:8,:)-foc([1;1],:);
p3=patch(xplot,yplot,[0.5 1 0.5],'DisplayName','Sextupoles');

ampl=0.03*tlim;     % Other multipoles
[xplot,yplot]=setxpl(sl(mpoles),ll(mpoles),[0;0;1;1;0],ampl*[0;1;1;0;0]);
%[xplot,yplot]=setxpl(sl(mpoles),ll(mpoles),[0;0;1;1;0],ampl*[-1;1;1;-1;-1]);
p4=patch(xplot,yplot,[0 0.5 0]);

ampl=0.015*tlim;    % BPMs
amplx=0.005*slim;
[xplot,yplot]=setxpl(sl(bpms),ones(1,sum(bpms)),[-amplx;0;amplx;0;-amplx],[0;-ampl;0;ampl;0]);
p5=patch(xplot,yplot,[0 0 0],'clipping','off');

if any(labs)
    slabs=sl(labs)+0.5*ll(labs);
    vlabs=cellfun(@(el) el.FamName,rok(labs),'UniformOutput',false);
    args={'Label',text(slabs,-0.03*tlim*ones(size(slabs)),vlabs,'Rotation',90,...
        'Interpreter','none','FontUnits','normalized','FontSize',0.025,...
        'HorizontalAlignment','right')};
else
    args={};
end
% Put patches in the background
set(ax,'Children',circshift(get(ax,'Children'),nlines));

pp=struct('Dipole',p1,'Quadrupole',p2,'Sextupole',p3,'Multipole',p4,...
    'BPM',p5,args{:});

    function [xpl,ypl]=setxpl(s,l,xmotif,ymotif)
        nm=length(xmotif);
        xpl=s(ones(nm,1),:)+l(ones(nm,1),:).*xmotif(:,ones(1,length(s)));
        ypl=ymotif(:,ones(1,length(s)));
    end

    function varargout=eltype(elem)
        if isfield(elem,'Class')
            cls=elem.Class;
        else
            cls=atguessclass(elem);
        end
        tst=strcmp(cls,{'Bend','Quadrupole','Sextupole','Multipole','Monitor'});
        varargout=num2cell(tst);
    end
end
