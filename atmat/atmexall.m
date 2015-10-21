function atmexall(varargin)
%ATMEXALL builds all AT platform deendent mex-files from C-sources
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required


PLATFORMOPTION = ['-D',computer,' ',sprintf('%s ',varargin{:})];
LIBDL='';
switch computer
case'GLNX86'
    LIBDL=' -ldl';
case 'GLNXA64'
    LIBDL=' -ldl';
end

try
    if ~verLessThan('matlab','7.11') %R2010b
        PLATFORMOPTION = [PLATFORMOPTION '-largeArrayDims '];
    end
    if ~verLessThan('matlab','8.3') %R2014a
        PLATFORMOPTION = [PLATFORMOPTION '-silent '];
    end
catch
end

% Navigate to the directory that contains pass-methods 
PASSMETHODDIR = fullfile(atroot,'..','atintegrators','');
mexpassmethod('all',PLATFORMOPTION);

% Navigate to the directory that contains tracking functions
cdir=fullfile(atroot,'attrack','');
MEXCOMMAND = ['mex ',PLATFORMOPTION,' -I''',PASSMETHODDIR,''' -outdir ',cdir,' ',fullfile(cdir,'atpass.c'),LIBDL];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% Navigate to the directory that contains some accelerator physics functions
cdir=fullfile(atroot,'atphysics','');
MEXCOMMAND = ['mex ',PLATFORMOPTION,' -outdir ',cdir,' -I''',PASSMETHODDIR,''' ',...
    fullfile(cdir,'findmpoleraddiffmatrix.c')];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% ADD 'MEXING' instructions for other C files
disp('ALL mex-files created successfully')
end
