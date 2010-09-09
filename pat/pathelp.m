% This provides usage information for the parallel version of AT.
% Xiaobiao Huang,  
% 7/12/2010
%
%% Introduction
%The existing mex-based AT passmethods are converted to having
%parallel-capability with the OpenMP approach. It enables the creation of
%multiple threads that splits the particle tracking work. The threads are
%run in parallel on different processors of a shared-memory computer, such
%as the usual multi-core PCs. This is useful for large tracking studies
%such as dynamic aperture and momentum aperture studies.

%% Compilation
%I have only tested the compilation on a Windows PC using VC 2008
%The commands to include OpenMP option looks like 
 mex DriftPass.c  CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" 
%The compilation flag may be different for different compilers.
%To compile all passmethods, you may use the script mexall.m
 mexall

%% Using the parallel passmethods

%To use the parallel capability, one needs to set the desired number of
%threads. This can be done by defining an environment variable
%OMP_NUM_THREADS. On my desktop (with 4 processors), I define the variable 
%to be 3. 
 getenv('OMP_NUM_THREADS')
 
 %You can set this variable dynamcially within matlab by the 'setenv'
 %command. 
 setenv('OMP_NUM_THREADS','4')
 getenv('OMP_NUM_THREADS')
 
 %% Test the parallel code
 
 sp3v82

 cd r:\
 which DriftPass
 X=zeros(6,1000); X(1,:) = 0.001*rand(1,1000);
tic; nX = ringpass(THERING, X, 100); toc    %old

 cd r:\xiahuang\misc\pat\test
 which DriftPass
 tic; pX = ringpass(THERING, X, 100); toc    %parallel

 norm(pX-nX)
 
 % The Windows task manager performance page shows the CPU work loads with
 % the sequential and parallel AT versions. For the test on my desktop, it
 % looks like 
 a=imread('task_mgr_omp.jpg');
 image(a)
