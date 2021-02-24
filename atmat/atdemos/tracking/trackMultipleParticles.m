function varargout=trackMultipleParticles(varargin)
%TRACKMULTIPLEPARTICLES     Time the tracking of multiple particles
%
%[TNORAD,TRAD]=TRACKMULTIPLEPARTICLES()
%   Track with default options, with and without longitudinal motion
%   and display the elapsed times.
%
%[...]=TRACKMULTIPLEPARTICLES(...,'nparticles',nparticles)
%   Set the number of particles (default 500)
%
%[...]=TRACKMULTIPLEPARTICLES(...,'nturns',nturns)
%   Set the number of turns (default 500)
%
%[...]=TRACKMULTIPLEPARTICLES(...,'omp_num_threads',nthreads)
%   Set the number of threads (Default: 0 (automatic) )

[nparticles,varargs]=getoption(varargin,'nparticles',500);
[nturns,varargs]=getoption(varargin,'nturns',500);
% Load the hmba lattice

a=load(fullfile(atroot,'../pyat/test_matlab/hmba.mat'));
mring=a.RING;
mring2=atradon(mring);

% Compute the initial parameters

dp=0.0;
[beamdata,~]=atx(mring, dp, 1);

% Build a x.x' particle distribution

rin=[atbeam(nparticles,beamdata.beam44(1:2,1:2));zeros(4,nparticles)];

% track without longitudinal motion
tic; rout1=ringpass(mring,rin,nturns,'silent',varargs{:}); t1=toc;
fprintf('%i particles, %i turns, without longitudinal motion: %f s\n', nparticles, nturns, t1);


% track with longitudinal motion

tic; rout2=ringpass(mring2,rin,nturns,'silent',varargs{:}); t2=toc;
fprintf('%i particles, %i turns,    with longitudinal motion: %f s\n', nparticles, nturns, t2);

if nargout == 0
    plot(rin(1,:),rin(2,:),'o'); hold on
    plot(rout1(1,:),rout1(2,:),'ro'); hold off
    xlabel('x [m]');
    ylabel('x'' [rd]');
    legend('Initial distribution', sprintf('After %i turns',nturns));
else
    varargout={t1, t2};
end