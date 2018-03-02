function h=plotContours(ring,indices,colorplotstr)
%PLOTCONTOURS Plots contours
%h=plotdriftcontours(ring)
%load S25ERRseed1.mat

[dat,prm]=atx(ring);
sig=cat(3,dat.beam66);

sigxx=squeeze(sig(1,1,:));
sigxy=squeeze(sig(1,3,:));
sigyy=squeeze(sig(3,3,:));

%dr278=findcells(ring,'FamName','DR_278');

hold on
for j=indices
    [x,y]=atmakeXYProjectionEllipse(sigxx(j),sigyy(j),sigxy(j));
    plot(x,y,colorplotstr)
end