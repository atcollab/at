%ATMEXALL builds all AT platform deendent mex-files from C-sources
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required


PLATFORMOPTION = ['-D',computer,' '];
LIBDL='';
switch computer
case'GLNX86'
    LIBDL=' -ldl';
case 'GLNXA64'
    LIBDL=' -ldl';
end

try
    if ~verLessThan('matlab','7.11')
        PLATFORMOPTION = [PLATFORMOPTION '-largeArrayDims '];
    end
catch
end

% Navigate to the directory that contains pass-methods 
cd(fullfile(atroot,'..','atintegrators',''));
PASSMETHODDIR = pwd;
disp(['Current directory: ',pwd]);
mexpassmethod('all',PLATFORMOPTION);

% Navigate to the directory that contains tracking functions
cd(fullfile(atroot,'attrack',''));
disp(['Current directory:', pwd]);
MEXCOMMAND = ['mex ',PLATFORMOPTION,'atpass.c',LIBDL];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% Navigate to the directory that contains some accelerator physics functions
cd(fullfile(atroot,'atphysics',''));
disp(['Current directory:', pwd]);
MEXCOMMAND = ['mex ',PLATFORMOPTION,' -I''',PASSMETHODDIR,''' findmpoleraddiffmatrix.c'];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% ADD 'MEXING' instructions for other C files
disp('ALL mex-files created successfully')
clear PASSMETHODDIR PLATFORMOPTION MEXCOMMAND
