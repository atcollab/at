function [pp,labels] = atplotsyn(ax,ring,options)
%ATPLOTSYN Helper function for ATPLOT
%
%PATCHES=ATPLOTSYN(AX,RING) Plots the magnetic elements found in RING

ylim=get(ax,'YLim');
slim=diff(get(ax,'XLim'));

sl=findspos(ring(:,1),1:size(ring,1)+1);
ll=diff(sl);

[dipoles,qpoles,spoles,mpoles,bpms]=cellfun(@eltype,ring);

ampl=0.05*ylim(2);
[xplot,yplot]=setxpl(sl(dipoles),ll(dipoles),[0;0;1;1;0],[0;ampl;ampl;0;0]);
p1=patch(xplot,yplot,[0.5 0.5 1]);

ampl=0.05*ylim(2);
foc=0.02*ylim(2)*sign(atgetfieldvalues(ring(qpoles),'PolynomB',{2}))';
[xplot,yplot]=setxpl(sl(qpoles),ll(qpoles),[0;0;0.5;1;1;0],[0;ampl;ampl;ampl;0;0]);
yplot(3,:)=yplot(3,:)+foc;
p2=patch(xplot,yplot,[1 0.5 0.5]);

ampl=0.04*ylim(2);
foc=0.01*ylim(2)*sign(atgetfieldvalues(ring(spoles),'PolynomB',{3}))';
[xplot,yplot]=setxpl(sl(spoles),ll(spoles),[0;0;0.33;0.66;1;1;0],[0;ampl;ampl;ampl;ampl;0;0]);
yplot(3:4,:)=yplot(3:4,:)+foc([1;1],:);
p3=patch(xplot,yplot,[0.5 1 0.5]);

ampl=0.03*ylim(2);
[xplot,yplot]=setxpl(sl(mpoles),ll(mpoles),[0;0;1;1;0],[0;ampl;ampl;0;0]);
p4=patch(xplot,yplot,[0 0.5 0]);

ampl=0.015*ylim(2);
amplx=0.005*slim;
[xplot,yplot]=setxpl(sl(bpms),ones(1,sum(bpms)),[-amplx;amplx;0;-amplx],[0;0;ampl;0]);
p5=patch(xplot,yplot,[0 0 0]);

pp=[p1;p2;p3;p4;p5];
labels={'Dipole','Quadrupole','Sextupole','Multipole','BPM'};

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
