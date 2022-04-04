function varargout=plClosedOrbit(varargin)
%PLCLOSEDORBIT Plots H and V 4x4 closed orbit
%
%Helper function for atplot: plot
%- H and V closed orbits on left axis
%- Dispersion on right axis
%
%  EXAMPLEs
% >> atbaseplot(ring,@plClosedOrbit,{'dp',0.01});
% >> atplot(ring,@plClosedOrbit,'dp',0.01);     (obsolete)
%
%  See also atplot atbaseplot

if nargout == 1 % From atplot
    lindata=varargin{1};
    CoD=cat(2,lindata.ClosedOrbit)';
    [xref,zref]=atreforbit(varargin{2});
    plotdata(1).values=CoD(:,[1 3])+[xref zref];    % left axis
    plotdata(1).labels={'x_{co}','z_{co}'};
    plotdata(1).axislabel='x,z [m]';
    dispersion=cat(2,lindata.Dispersion)';
    plotdata(2).values=dispersion(:,3);             % right axis
    plotdata(2).labels={'\eta_y'};
    plotdata(2).axislabel='Vertical dispersion [m]';
    varargout={plotdata};
else % From atbaseplot
    [ring,dpp]=deal(varargin{1:2});
    [linargs,varargs]=opticsoptions(varargin(3:end));
    [~,lindata]=atlinopt6(ring,1:length(ring)+1,linargs{:});
    s=cat(1,lindata.SPos);
    varargout={s,plClosedOrbit(lindata,ring,dpp,varargs{:})};
end
end
