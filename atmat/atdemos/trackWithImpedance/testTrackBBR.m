%create a fast ring from ESRF lattice
ring=atradoff(atsetcavity(esrf(),8e6,1,992));

%Now generate fast ring and add impedance element
[fastring,fastringrad]=atfastring(ring);

%generate a Broad Band Resonator Impedance
impedance_elem = makeBBR(findspos(ring, length(ring)+1) , 0.005);

%add to fast ring
fastringBBR=[fastringrad;impedance_elem];

%Generate particle distribution and track
%watch rms values
NP=1e4;
Nturns=1000;
ring=atradoff(esrf());
particles0=atbeam(NP,atsigma(ring));
figure(1);
% histogram(particles0(6,:));
hist(particles0(6,:), 35);
xlabel('s [m]');
ylabel('# particles');
title('Initial distribution');

particles=ringpass(fastringBBR, particles0, Nturns, 'Silent');

figure(2);
%histogram(particles(6,:));
hist(particles(6,:),35);
xlabel('s [m]');
ylabel('# particles');
title('After 1000 turns');
