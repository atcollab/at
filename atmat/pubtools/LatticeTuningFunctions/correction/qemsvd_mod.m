function dq=qemsvd_mod(a,b,neig,plot)
%function dq=qemsvd_mod(a,b,neig,plot)
% given response matrix a and vector b to be corrected the function finds 
% the vector dq so that b-a*dq=0
%
% if the input plot is = 1 the plot of correction effect and correctors
% strengths versus number of eigenvectors is also computed.
%
% svd part is carbon copy of qemsvd by L.Farvacque in qempanel

[u,s,v]=svd(a,0);
lambda=diag(s);
nmax=length(lambda);
eigsorb=u'*b;
if neig > nmax
    neig=nmax;
    warning('Svd:maxsize',['number of vectors limited to ' num2str(nmax)]);
end
eigscor=eigsorb(1:neig)./lambda(1:neig);
dq=v(:,1:neig)*eigscor;


% plot correction effect to set appropriate number of eigenvectors.
if nargin<4
    plot=0;
end

if plot
    numeigen0=neig;
    neig=sort([1:5:nmax,numeigen0]);
    dqstd=zeros(size(neig));
    ostd=zeros(size(neig));
    for ineig=1:length(neig)% loop number of eigenvectors used in the correction
        eigscor=eigsorb(1:neig(ineig))./lambda(1:neig(ineig));
        dqe=v(:,1:neig(ineig))*eigscor;
        o=b-a*dqe;
        dqstd(ineig)=std(dqe); 
        ostd(ineig)=std(o);
    end
    
    f=figure('visible','on');
    AX=plotyy(neig,ostd,neig,dqstd);
    xlabel('number of eigenvectors');
    set(get(AX(1),'Ylabel'),'String','rms vector to correct')
    set(get(AX(2),'Ylabel'),'String','rms correctors')
    %set(AX(1),'YScale','log');
    set(AX(2),'YScale','log');
    title('Correctors strength and residual vs number of eigenvectors')
    saveas(gca,[datestr(now,'yyyymmddTHHMMSS') '_eigplot.fig']);
    export_fig([datestr(now,'yyyymmddTHHMMSS') '_eigplot.jpg']);
    close(f);
    
end
