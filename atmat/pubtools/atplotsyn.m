function pp = atplotsyn(ax,ring)
%ATPLOTSYN Helper function for ATPLOT
%
%PATCHES=ATPLOTSYN(AX,RING) Plots the magnetic elements found in RING

ylim=get(ax,'YLim');

magnets=getelems(ring,'PolynomB');
sl=findspos(ring,1:length(ring)+1);
ll=diff(sl);
b3=zeros(1,length(ring));
b3(magnets)=getcellstruct(ring,'PolynomB',find(magnets),3);
b2=zeros(1,length(ring));
b2(magnets)=getcellstruct(ring,'PolynomB',find(magnets),2);

dipoles=getelems(ring,'BendingAngle');
dips=cat(1,ring{dipoles});
l=ll(dipoles);
s=sl(dipoles);
lplot=zeros(5,length(l));
lplot(3:4,:)=[l;l];
xplot=s(ones(5,1),:)+lplot;
yplot=zeros(5,length(l));
yplot(2:3,:)=0.05*ylim(2);
%p1=patch(xplot,yplot,[0.5 0.5 1],'FaceAlpha',0.2);
p1=patch(xplot,yplot,[0.5 0.5 1]);

qpoles=(b2~=0) & ~dipoles;
strength=b2(qpoles);
l=ll(qpoles);
s=sl(qpoles);
lplot=zeros(6,length(l));
lplot(3:5,:)=[0.5*l;l;l];
xplot=s(ones(6,1),:)+lplot;
yplot=zeros(6,length(l));
yplot(2:4,:)=0.05*ylim(2);
yplot(3,:)=yplot(3,:)+0.02*ylim(2)*sign(strength);
p2=patch(xplot,yplot,[1 0.5 0.5]);

spoles=(b3~=0) & ~dipoles & ~qpoles;
strength=b3(spoles);
l=ll(spoles);
s=sl(spoles);
lplot=zeros(6,length(l));
lplot(3:5,:)=[0.5*l;l;l];
xplot=s(ones(6,1),:)+lplot;
yplot=zeros(6,length(l));
yplot(2:4,:)=0.03*ylim(2);
yplot(3,:)=yplot(3,:)+0.01*ylim(2)*sign(strength);
p3=patch(xplot,yplot,[0.5 1 0.5]);
pp=[p1;p2;p3];

    function res=getelems(cellarray,fieldname)
        res=false(1,length(cellarray));
        for i=1:length(cellarray)
            res(i)=isfield(cellarray{i},fieldname);
        end
    end
        
end

