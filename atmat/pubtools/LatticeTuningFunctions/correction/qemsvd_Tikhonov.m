function dq=qemsvd_Tikhonov(a,b,lambda,plot)
%function dq=qemsvd_Tikhonov(a,b,lambda)
% given response matrix a and vector b to be corrected the function finds 
% the vector dq so that b-a*dq=0
% lambda is the tikhonov parameter
%
%see also: svd, tikhonov

[u,s,v]=csvd(a,0);
x_0=zeros(size(b)); % initial correction vector
[dq,~,~] = tikhonov(u,s,v,b,lambda);

% plot correction effect to set appropriate number of eigenvectors.
if nargin<4
    plot=0;
end

if plot
    lambda0=lambda;
    neig=sort([0.1:0.1:lambda0,lambda0+0.1:0.25:3]);
    dqstd=zeros(size(neig));
    ostd=zeros(size(neig));
    for ineig=1:length(neig)% loop number of eigenvectors used in the correction
        [dqe,~,~] = tikhonov(u,s,v,b,neig(ineig));
        o=b-a*dqe;
        dqstd(ineig)=std(dqe); 
        ostd(ineig)=std(o);
    end
    
    f=figure('visible','on');
    AX=plotyy(neig,ostd,neig,dqstd);
    xlabel('Tikhonov parameter');
    set(get(AX(1),'Ylabel'),'String','rms vector to correct')
    set(get(AX(2),'Ylabel'),'String','rms correctors')
    %set(AX(1),'YScale','log');
    set(AX(2),'YScale','log');
    title('Correctors strength and residual vs Tikhonov parameter')
    saveas(gca,[datestr(now,'yyyymmddTHHMMSS') '_Tikhonovplot.fig']);
    export_fig([datestr(now,'yyyymmddTHHMMSS') '_Tikhonovplot.jpg']);
    close(f);
    
end
end

