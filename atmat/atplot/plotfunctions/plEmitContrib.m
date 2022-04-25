function varargout=plEmitContrib(varargin)
%PLEMITCONTRIB  Plot H/rho³ at every dipole
%
% USAGE:
% >> atbaseplot(ring,@PLEMITCONTRIB);
% >> atplot(ring,@PLEMITCONTRIB);     (obsolete)
%
%  See also atbaseplot

if nargout == 1 % From atplot
    [lindata,ring,dpp]=deal(varargin{1:3});
    idx=1:length(ring)+1;
    H=CurlyH(ring,dpp,idx(:)');

    r=zeros(size(H));
    bend=findcells(ring,'BendingAngle');
    r(bend)=getcellstruct(ring,'Length',bend)./getcellstruct(ring,'BendingAngle',bend);
    emitcontr=H./(r.^3)*1e9;
    emitcontr(isinf(emitcontr))=0;

    beta=cat(1,lindata.beta);                     % left axis
    plotdata(1).values=beta;
    plotdata(1).labels={'\beta_x','\beta_z'};
    plotdata(1).axislabel='\beta [m]';

    dispersion=cat(2,lindata.Dispersion)'; % right axis
    plotdata(2).values=[dispersion(:,1)*100 emitcontr H*10000];
    plotdata(2).labels={'\eta_x cm','H/r^{3}*1e9','H 10-4'};
    plotdata(2).axislabel='dispersion [cm]';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=opticsoptions(varargin(3:end));
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    s=cat(1,lindata.SPos);
    varargout={s,plEmitContrib(lindata,ring,dpp,varargs{:})};
end

end
