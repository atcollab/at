%give some initial coordinates.  Track through sample lattices.
%Store results for tracking tests.
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

Z1_dba=ringpass(dba_ring,Z0,1);
Z1_FODO=ringpass(FODO_ring,Z0,1);
Z1_AS=ringpass(AS_ring,Z0,1);
Z1_esrf=ringpass(esrf_ring,Z0,1);
Z1_soleil=ringpass(soleil_ring,Z0,1);
Z1_thomx=ringpass(thomx_ring,Z0,1);

save('trackTestData.mat','Z0','Z1_dba','Z1_FODO','Z1_AS','Z1_esrf','Z1_soleil','Z1_thomx')