function plotdata=plSigmaSigmap(lindata,ring,dpp,varargin) 
%PLSIGMASIGMAP Plots beam sizes and divergences
% Plots sigmax and sigmay on left axis and 
%       sigmax' and sigmay' on right axis

pos=1:length(lindata);

[lin,pm]=atx(ring,dpp,pos); % @ ID (element 1)

ex=arrayfun(@(x)x.emit44(1),lin);
ey=arrayfun(@(x)x.emit44(2),lin);

% dx@id
dx=arrayfun(@(x)x.Dispersion(1),lindata);
dxp=arrayfun(@(x)x.Dispersion(2),lindata);
dy=arrayfun(@(x)x.Dispersion(3),lindata);
dyp=arrayfun(@(x)x.Dispersion(4),lindata);
bx=arrayfun(@(x)x.beta(1),lindata);
ax=arrayfun(@(x)x.alpha(1),lindata);
by=arrayfun(@(x)x.beta(2),lindata);
ay=arrayfun(@(x)x.alpha(2),lindata);

% sigmaX*sigmaXp @ id

sx=sqrt(bx.*ex+(pm.espread.*dx).^2);
sxp=sqrt((1+ax.^2)./bx.*ex+(pm.espread.*dxp).^2);
%sxsxp=sx.*sxp;

sy=sqrt(by.*ey+(pm.espread.*dy).^2);
syp=sqrt((1+ay.^2)./by.*ey+(pm.espread.*dyp).^2);
%sysxp=sy.*syp;

plotdata(1).values=[sx; sy]';
plotdata(1).labels={'\sigma_{x}','\sigma_{y}'};
plotdata(1).axislabel=' and  [m]';

plotdata(2).values=[sxp; syp]';
plotdata(2).labels={'\sigma_{x}''','\sigma_{y}'''};
plotdata(2).axislabel=' [m]';

end