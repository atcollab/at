function varargout=plotbetadisp(varargin)
%function [s,plotdata]=plotbetadisp(ring,dpp,plotfun,varargin)
%PLOTBETADISP Plot beta functions and dispersion
%
%Helper function for tplot:
%- beta functions on the left axis
%- diespersion on the right axis
%
%  USAGE:
% >> atbaseplot(ring,@plotbetadisp);
% >> atplot(ring,@plotbetadisp);      (obsolete)
%
%See also atplot atbaseplot

if nargout == 1     % From atplot
    lindata=varargin{1};
    beta=cat(1,lindata.beta);
    plotdata(1).values=beta;
    plotdata(1).labels={'\beta_x','\beta_z'};
    plotdata(1).axislabel='\beta [m]';
    dispersion=cat(2,lindata.Dispersion)';
    plotdata(2).values=dispersion(:,1);
    plotdata(2).labels={'\eta_x'};
    plotdata(2).axislabel='dispersion [m]';
    varargout={plotdata};
else                % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=linoptions(getdparg(varargin(2:end)));
    [ringdata,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:}); %#ok<ASGLU>
    s=cat(1,lindata.SPos);
    varargout={s,plotbetadisp(lindata,ring,dpp,varargs{:})};
end
end
