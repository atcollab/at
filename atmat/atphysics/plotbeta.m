function varargout = plotbeta(varargin)
%PLOTBETA plots UNCOUPLED! beta-functions
% PLOTBETA(RING) calculates beta functions of the lattice RING
% PLOTBETA with no argumnts uses THERING as the default lattice
%  Note: PLOTBETA uses FINDORBIT4 and LINOPT which assume a lattice
%  with NO accelerating cavities and NO radiation
%
% See also PLOTCOD 
global THERING
if nargin == 0
	RING = THERING;
else
    RING = varargin{1};
end

L = length(RING);
spos = findspos(RING,1:L+1);

[TD, tune] = twissring(RING,0,1:(length(RING)+1));
BETA = cat(1,TD.beta);
S  = cat(1,TD.SPos);

disp(tune)

if nargin > 1 & ishandle(varargin{2})
    figure(varargin{2});
else
    figure
end
% plot betax and betay in two subplots

subplot(2,1,1)
plot(S,BETA(:,1),'.-b');

A = axis;
A(1) = 0;
A(2) = S(end);
axis(A);
%xlabel('s - position [m]');
ylabel('\beta_x [m]');
grid on


title('\beta-functions');

subplot(2,1,2)
plot(S,BETA(:,2),'.-r');
% Set the same horizontal scale as beta x plot
B = axis;
axis([A(1:2) B(3:4)]);
xlabel('s - position [m]');
ylabel('\beta_y [m]');
grid on
