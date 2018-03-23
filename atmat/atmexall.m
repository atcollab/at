function atmexall(varargin)
%ATMEXALL builds all AT platform dependent mex-files from C-sources
%
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required

% Add AT directories to the Matlab path
atpath

PLATFORMOPTION = ['-D',computer,' ',sprintf('%s ',varargin{:})];
LIBDL='';
switch computer
    case'GLNX86'
        LIBDL=' -ldl';
    case 'GLNXA64'
        LIBDL=' -ldl';
end

try
    if ~verLessThan('matlab','9.4')
        PLATFORMOPTION = [PLATFORMOPTION '-R2018a '];
    elseif ~verLessThan('matlab','7.11') %R2010b
        PLATFORMOPTION = [PLATFORMOPTION '-largeArrayDims '];
    end
    if ~verLessThan('matlab','8.3') %R2014a
        PLATFORMOPTION = [PLATFORMOPTION ' -silent '];
    end
catch
end

% Navigate to the directory that contains tracking functions
lastwarn('');

PASSMETHODDIR = fullfile(atroot,'..','atintegrators','');
cdir=fullfile(atroot,'attrack','');
MEXCOMMAND = ['mex ',PLATFORMOPTION,'-outdir ',cdir,' -I''',PASSMETHODDIR,''' ',fullfile(cdir,'atpass.c'),LIBDL];
disp(MEXCOMMAND);
eval(MEXCOMMAND);
[warnmess,warnid]=lastwarn; %#ok<ASGLU>
if strcmp(warnid,'MATLAB:mex:GccVersion_link')
    warning('Disabling the compiler warning');
end

% Navigate to the directory that contains some accelerator physics functions
oldwarns=warning('OFF','MATLAB:mex:GccVersion_link');
cdir=fullfile(atroot,'atphysics','Radiation');
MEXCOMMAND = ['mex ',PLATFORMOPTION,'-outdir ',cdir,' -I''',PASSMETHODDIR,''' ',fullfile(cdir,'findmpoleraddiffmatrix.c')];
disp(MEXCOMMAND);
eval(MEXCOMMAND);
cdir=fullfile(atroot,'atphysics','NonLinearDynamics');
MEXCOMMAND = ['mex ',PLATFORMOPTION,'-outdir ',cdir,' ',fullfile(cdir,'RDTelegantAT.cpp')];
disp(MEXCOMMAND);
eval(MEXCOMMAND);

% NAFF
cdir=fullfile(atroot,'atphysics','nafflib');
MEXCOMMAND = ['mex ',PLATFORMOPTION,'-outdir ',cdir,' ',fullfile(cdir,'nafflib.c'),' ',...
    fullfile(cdir,'modnaff.c'),' ',...
    fullfile(cdir,'complexe.c'),' ',...
    ];
disp(MEXCOMMAND);
eval(MEXCOMMAND);


% Navigate to the directory that contains pass-methods 
mexpassmethod('all',PLATFORMOPTION);
warning(oldwarns.state,oldwarns.identifier);

% ADD 'MEXING' instructions for other C files
%disp('ALL mex-files created successfully')
end
