function varargout=plBeamSize(varargin)
%PLBEAMSIZE Plot H and V beam size
%
% USAGE:
% >> atbaseplot(ring,@PLBEAMSIZE);
% >> atplot(ring,@PLBEAMSIZE);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [lindata,ring,dpp]=deal(varargin{1:3}); %#ok<ASGLU> 
    [r2,radindex,~]=atradon(ring);
    [envelope,~,~]=ohmienvelope(r2,radindex,1:length(r2)+1);

    plotdata(1).values=cat(1,lindata.beta);         % left axis
    plotdata(1).labels={'\beta_x','\beta_z'};
    plotdata(1).axislabel='\beta [m]';

    sigma=real(cat(2,envelope.Sigma))*1e3;
    dispersion=cat(2,lindata.Dispersion)';

    plotdata(2).values=[dispersion(:,1) sigma'];    % right axis
    plotdata(2).labels={'\eta_x ','\sigma_x [mm]','\sigma_y [mm]'};
    plotdata(2).axislabel='dispersion [m] and \sigma [mm]';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    s=findspos(ring,1:length(ring)+1);
    [linargs,vargs] = linoptions(varargin(4:end),dpp);
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    varargout={s,plBeamSize(lindata,ring,dpp,vargs{:})};
end
end
