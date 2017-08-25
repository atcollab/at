function varargout = plotcod(RING,DP)
%PLOTCOD Closed Orbit Distortion
% PLOTCOD(RING,DP) finds and plots closed orbit for a given momentum 
%  deviation DP. It calls FINDORBIT4 which assumes a lattice
%  with NO accelerating cavities and NO radiation

localspos = findspos(RING,1:length(RING)+1);
orbit = findorbit4(RING,DP,1:length(RING)+1);


plot(localspos,orbit(1,:),'.-r');
title('Closed Orbit Distortion')
hold on
plot(localspos,orbit(3,:),'.-b');
hold off

A = axis;
A(1) = 0;
A(2) = localspos(end);
axis(A);

legend('Horizontal','Vertical');
xlabel('s - position [m]');
ylabel('orbit [m]');

grid on

if nargout > 0
	varargout{1} = orbit;
end