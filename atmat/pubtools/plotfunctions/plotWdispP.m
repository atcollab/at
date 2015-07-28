function varargout=plotWdispP(varargin)
%plotWdispP    Plot W functions
%
%Helper function for atplot: plot
%- W functions (derivatives of beta-functions versus momentum) on left axis
%- derivative of dispersion on right axis

if nargout == 1
    [ring,dpp]=deal(varargin{2:3});
    refpts=1:length(ring)+1;
    DE=0.001;
    
    [lz,~,~]=atlinopt(ring,dpp,refpts);
    [lpd,~,~]=atlinopt(ring,dpp+DE,refpts);
    [lmd,~,~]=atlinopt(ring,dpp-DE,refpts);
    bz=cat(1,lz.beta);
    bp=cat(1,lpd.beta);
    bm=cat(1,lmd.beta);
    az=cat(1,lz.alpha);
    ap=cat(1,lpd.alpha);
    am=cat(1,lmd.alpha);
    
    aa=((bp-bm)./bz)/2/DE;
    bb=((ap-am)-az./bz.*(bp-bm))/2/DE;
    
    plotdata(1).values=sqrt(aa.^2+bb.^2);
    plotdata(1).labels={'W_x','W_z'};%{'\beta_x/D\delta','\beta_z/D\delta'};
    plotdata(1).axislabel='W_x W_z [m]';
    
    dp=cat(2,lpd.Dispersion)';
    dm=cat(2,lmd.Dispersion)';
    plotdata(2).values=100*(dp(:,1)-dm(:,1))/2/DE;
    plotdata(2).labels={'\partial\eta_{x}/\partial\delta'};
    plotdata(2).axislabel='D'' [cm]';
    varargout={plotdata};
else
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plotWdispP([],varargin{:})};
end
end
