@echo off
echo.
echo Starting MATLAB and running Frequency Map Analysis...
echo Input file : 
echo        %%CD%%aspsr_v1_fma_input.dat
echo Matlab error logs saved to : 
echo        %%CD%%aspsr_v1_fma_input_matlaberrors.log
echo.
matlab -nosplash -nodesktop -nojvm -logfile %%CD%%aspsr_v1_fma_input_matlaberrors.log -r "fma = fma_lib('analysis','%%CD%%aspsr_v1_fma_input.dat'); plot_emailresults(fma); exit;"
