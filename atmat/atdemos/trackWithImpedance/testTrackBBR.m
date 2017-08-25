%create a fast ring from ESRF lattice
%generate a Broad Band Resonator Impedance
%add to fast ring
makeFastRingWithBBR;

%Generate particle distribution and track
%watch rms values
NP=1e4;
Nturns=1000;
particles=atbeam(NP,atsigma(esrf));
disp(std(particles'))
for i=1:Nturns
  particles=ringpass(fastringBBR,particles,1);
  disp(std(particles'))
end