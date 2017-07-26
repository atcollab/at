% Demo of using atfastring for fast tracking
% The ESRF lattice is loaded, cavity is set, and then the fast ring is
% created.
% An electron is tracked through the full lattice and fast ring and
% the tracking time and tracking results are compared. 

clear all
clc
echo on

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


% Check tunes and chromaticity:
[p,t,c]=atlinopt(ring,0,1);
[pf,tf,cf]=atlinopt(fastring,0,1);

t
tf
c
cf

% Now compare the horizontal and transverse tracking results.
plot(1:length(z1),z1(1,:),'r',1:length(z1fast),z1fast(1,:),'b');
legend('full ring','fast ring');
xlabel('turns');
ylabel('x (m)');

figure
plot(z1(5,:),'r');
hold on
plot(z1fast(5,:),'b');
legend('full ring','fast ring');
xlabel('turns');
ylabel('delta');


% figure
% hold on
% plot(z1fastrad(5,:),'-k');