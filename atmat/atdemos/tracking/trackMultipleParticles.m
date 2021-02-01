function trackMultipleParticles(varargin)
%TRACKMULTIPLEPARTICLES     Time the tracking of multiple particles
%
%TRACKMULTIPLEPARTICLES()   Track with default options, with and without
%                           longitudinal motion and display the elapsed time
%
%TRACKMULTIPLEPARTICLES(...,'nparticles',nparticles)
%   Set the number of particles (default 500)
%
%TRACKMULTIPLEPARTICLES(...,'nturns',nturns)
%   Set the number of turns (default 500)
     
[nparticles,varargs]=getoption(varargin,'nparticles',500);
[nturns,~]=getoption(varargin,'nturns',500);
% Load the hmba lattice

a=load(fullfile(atroot,'../pyat/test_matlab/hmba.mat'));
mring=a.RING;
mring2=atradon(mring);

% Compute the initial parameters

dp=0.0;
[beamdata,~]=atx(mring, dp, 1);

% Build a x.x' particle distribution

rin=[atbeam(nparticles,beamdata.beam44(1:2,1:2));zeros(4,nparticles)];
figure(1);
plot(rin(1,:),rin(2,:),'o'); hold on

% track without longitudinal motion
tic; rout=ringpass(mring,rin,nturns,'silent'); t=toc;
fprintf('%i particles, %i turns, without longitudinal motion: %f s\n', nparticles, nturns, t);

plot(rout(1,:),rout(2,:),'ro'); hold off
xlabel('x [m]');
ylabel('x'' [rd]');
legend('Initial distribution', sprintf('After %i turns',nturns));

% track with longitudinal motion

tic; rout=ringpass(mring2,rin,nturns,'silent'); t=toc;
fprintf('%i particles, %i turns,    with longitudinal motion: %f s\n', nparticles, nturns, t);

end