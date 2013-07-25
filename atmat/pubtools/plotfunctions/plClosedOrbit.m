function plotdata=plClosedOrbit(lindata,ring,dpp,varargin) %#ok<INUSD>
%DEFAULTPLOT    Default plotting function for ATPLOT
%
%Plots beta-functions on left axis and dispersion on right axis

CoD=cat(2,lindata.ClosedOrbit);
plotdata(1).values=[CoD(1,:)' CoD(3,:)'];
plotdata(1).labels={'x','z'};
plotdata(1).axislabel='y_{co} and x_{co} [m]';
dispersion=cat(2,lindata.Dispersion)';
plotdata(2).values=dispersion(:,3);
plotdata(2).labels={'\eta_y'};
plotdata(2).axislabel='vertical dispersion [m]';

end
