function plotdata=pltmisalignments(lindata,ring,dpp,varargin) %#ok<INUSD>
%DEFAULTPLOT    Default plotting function for ATPLOT
%
%Plots orbits, misalignments and tilts.

CoD=cat(2,lindata.ClosedOrbit);
[X,Y,~,Tilt]=GetExistingErrors(ring);
% T=T-pi/2;
% T(T==pi/2)=0;
% T(T==-pi/2)=0;

X=[X,0];
Y=[Y,0];
Tilt=[Tilt,0];
% tiltedelem=findcells(ring,'Tilt');
% tiltval=[getcellstruct(ring,'Tilt',tiltedelem); 0];

% plotdata(1).values=[CoD(1,:)' CoD(3,:)' X' Y']*1e6;%
% plotdata(1).labels={'x','z','x misal','z misal'};
% plotdata(1).axislabel='orbits and misalignments [micro m]';
plotdata(1).values=[X' Y']*1e6;%
plotdata(1).labels={'x misal','z misal'};
plotdata(1).axislabel='misalignments [micro m]';
% % CoD(1,:)' 
plotdata(2).values=[Tilt']*1e6;
plotdata(2).labels={'Rot. s-axis'};
plotdata(2).axislabel='Rotation [micro rad]';

end 