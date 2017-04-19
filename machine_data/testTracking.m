function out = testTracking
%test results of tracking against values previously computed
%trackTestData.mat should already exist.  This tests to see
%if tracking results are unchanged.
xextent = 0.001;
yextent = 0.001;
n = 10;
Z0=[];
for j=-n:n
    for k=-n:n
        x0=j*xextent/n;
        y0=j*yextent/n;
        PS0=[x0;0;y0;0;0;0];
        Z0=[Z0,PS0];
    end
end
%Now, load each lattice and track
dba_ring=dba;
FODO_ring=FODO;
AS_ring=australian_synchrotron;
esrf_ring=esrf;
soleil_ring=soleil;
thomx_ring=thomx;

Z1_dba_=ringpass(dba_ring,Z0,1);
Z1_FODO_=ringpass(FODO_ring,Z0,1);
Z1_AS_=ringpass(AS_ring,Z0,1);
Z1_esrf_=ringpass(esrf_ring,Z0,1);
Z1_soleil_=ringpass(soleil_ring,Z0,1);
Z1_thomx_=ringpass(thomx_ring,Z0,1);

load trackTestData.mat

epsilon=1e-9;
% out = [isequal(Z1_dba_,Z1_dba),isequal(Z1_FODO_,Z1_FODO),isequal(Z1_AS_,Z1_AS),...
%     isequal(Z1_esrf_,Z1_esrf),isequal(Z1_soleil_,Z1_soleil),isequal(Z1_thomx_,Z1_thomx)];
out=[max(max(Z1_dba_-Z1_dba))<epsilon,max(max(Z1_FODO_-Z1_FODO))<epsilon,...
    max(max(Z1_AS_-Z1_AS))<epsilon,max(max(Z1_esrf_-Z1_esrf))<epsilon,...
    max(max(Z1_soleil_-Z1_soleil))<epsilon,max(max(Z1_thomx_-Z1_thomx))<epsilon];