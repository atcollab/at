%Demo of using atfastring for fast tracking
%load esrf ring
ring=esrf;

% indcav=findcells(ring,'Class','RFCavity');
% cav=ring(indcav(1));
% ring(indcav(:))=[];
% ring=[cav;ring];

% set RF cavity
ring=atsetcavity(ring,8e6,0,992);

%Now, create fastring and fastringrad from ring.
[fastring,fastringrad]=atfastring(ring);


%Set an initial condition for tracking
z0=[1e-5,0,0,0,1e-3,0]';

%Now, track with full ring and with fast ring and time the computation.
tic
z1=ringpass(ring,z0,500);
toc

tic
z1fast=ringpass(fastring,z0,500);
toc


% Now compare the horizontal and transverse tracking results.
plot(z1(1,:),'r');
hold on
plot(z1fast(1,:),'b');

figure
plot(z1(5,:),'r');
hold on
plot(z1fast(5,:),'b');

% figure
% hold on
% plot(z1fastrad(5,:),'-k');