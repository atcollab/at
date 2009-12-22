%ATMEXALL builds all AT platform deendent mex-files from C-sources
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required



ATROOT = atroot
% Navigate to the directory that contains pass-methods 
cd(ATROOT)
cd simulator
cd element
PASSMETHODDIR = pwd;
disp(['Current directory: ',PASSMETHODDIR]);
mexpassmethod('all');

% Navigate to the directory that contains tracking functions
cd(ATROOT)
cd simulator
cd track

disp(['Current directory:', pwd]);


PLATFORMOPTION = ['-D',computer,' '];

MEXCOMMAND = ['mex ',PLATFORMOPTION,'atpass.c'];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% Navigate to the directory that contains some accelerator physics functions
cd(ATROOT)
cd atphysics
disp(['Current directory:', pwd]);

% findmpoleraddiffmatrix.c
disp('mex findmpoleraddiffmatrix.c')

eval(['mex findmpoleraddiffmatrix.c -I''',PASSMETHODDIR,'''']);
% ADD 'MEXING' instructions for other C files
disp('ALL mex-files created successfully')
clear ATROOT PASSMETHODDIR WARNMSG PLATFORMOPTION MEXCOMMAND
