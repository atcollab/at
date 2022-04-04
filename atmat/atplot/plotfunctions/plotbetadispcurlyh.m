function varargout=plotbetadispcurlyh(varargin)
%PLOTBETADISPCURLYH Plot beta, dispersion and H
%
% USAGE:
% >> atbaseplot(ring,@PLOTBETADISPCURLYH);
% >> atplot(ring,@PLOTBETADISPCURLYH);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [lindata,~,~]=deal(varargin{1:3});
    H=CurlyHlindata(lindata);

    beta=cat(1,lindata.beta);                     % left axis
    plotdata(1).values=beta;
    plotdata(1).labels={'\beta_x','\beta_z'};
    plotdata(1).axislabel='\beta [m]';
    dispersion=cat(2,lindata.Dispersion)'; % right axis

    plotdata(2).values=[dispersion(:,1)*100 H*10000];
    plotdata(2).labels={'\eta_x','H*10^{2}'};
    plotdata(2).axislabel='\eta_x, H [cm] ';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=opticsoptions(varargin(3:end));
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    s=cat(1,lindata.SPos);
    varargout={s,plotbetadispcurlyh(lindata,ring,dpp,varargs{:})};
end
end
