function plotdata=plenvelope(lindata,ring,dpp,varargin) %#ok<INUSD>
elt=length(ring);
lind=atx(ring,dpp,1:elt+1);
sigma=cat(3,lind.beam66);
plotdata(1).values=sqrt([squeeze(sigma(1,1,:)) squeeze(sigma(3,3,:))]);
plotdata(1).labels={'horizontal','vertical'};
plotdata(1).axislabel='beam envelope';
end
