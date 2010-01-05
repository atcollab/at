%ATMEXALL builds all AT platform deendent mex-files from C-sources
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required


PLATFORMOPTION = ['-D',computer,' '];
LIBDL='';
switch computer
case 'GLNX86'
    PLATFORMOPTION = [PLATFORMOPTION,'CC=gcc4 LD=gcc4 ']; 
    LIBDL=' -ldl';
case 'GLNXA64'
    LIBDL=' -ldl';
end

ATROOT = atroot
% Navigate to the directory that contains pass-methods 
cd(ATROOT)
cd simulator
cd element
PASSMETHODDIR = pwd;
disp(['Current directory: ',PASSMETHODDIR]);
mexpassmethod('all',PLATFORMOPTION);

% Navigate to the directory that contains tracking functions
cd(ATROOT)
cd simulator
cd track

disp(['Current directory:', pwd]);

MEXCOMMAND = ['mex ',PLATFORMOPTION,'atpass.c',LIBDL];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% Navigate to the directory that contains some accelerator physics functions
cd(ATROOT)
cd atphysics
disp(['Current directory:', pwd]);

% findmpoleraddiffmatrix.c
MEXCOMMAND = ['mex ',PLATFORMOPTION,'findmpoleraddiffmatrix.c -I''',PASSMETHODDIR,''''];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% ADD 'MEXING' instructions for other C files
disp('ALL mex-files created successfully')
clear ATROOT PASSMETHODDIR WARNMSG PLATFORMOPTION MEXCOMMAND
