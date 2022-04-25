function varargout=plotsqrtbetadispcurlyh(varargin)
%PLOTSQRTBETADISPCURLYH Plot sqrt(beta), dispersion and H
%
% USAGE:
% >> atbaseplot(ring,@PLOTSQRTBETADISPCURLYH);
% >> atplot(ring,@PLOTSQRTBETADISPCURLYH);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [lindata,~,~]=deal(varargin{1:3});
    H=CurlyHlindata(lindata);

    beta=cat(1,lindata.beta);                     % left axis
    plotdata(1).values=sqrt(beta);
    plotdata(1).labels={'\beta_x^{0.5}','\beta_z^{0.5}'};
    plotdata(1).axislabel='\beta^{0.5} [m^{0.5}]';
    dispersion=cat(2,lindata.Dispersion)'; % right axis

    plotdata(2).values=[dispersion(:,1)*100 H*10000];
    plotdata(2).labels={'\eta_x [cm]','H [10-4]'};
    plotdata(2).axislabel='dispersion [cm]';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=opticsoptions(varargin(3:end));
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    s=cat(1,lindata.SPos);
    varargout={s,plotsqrtbetadispcurlyh(lindata,ring,dpp,varargs{:})};
end
end
