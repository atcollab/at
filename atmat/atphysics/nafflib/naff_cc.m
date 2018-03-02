function naff_cc
%NAFF_CC Compile nafflibrary for Matlab
%

% Modified by Laurent S. Nadolski
% April 6th, 2007

cd_old = pwd;
cd(fileparts(which('naff_cc')))

disp(['Compiling NAFF routines on ', computer,'.'])

switch computer
    case 'SOL2'
        PLATFORMOPTION = ['-D',computer,' '];
    case 'GLNXA64'
        PLATFORMOPTION = ['-ldl -D',computer,' ']; % added by Laurent April 6th, 2007        
    case 'GLNX86' 
        PLATFORMOPTION = ['-ldl -D',computer,' ']; % added by Laurent April 6th, 2007
    case  'PCWIN64'
        PLATFORMOPTION = ['-D',computer,' LDFLAGS=''-pthread -shared -m64'' '];
    case  'PCWIN'
        PLATFORMOPTION = ['-D',computer,' '];
    case  'MACI64'
        PLATFORMOPTION = ['-D',computer,' LDFLAGS=''-pthread -shared -m64'' '];
    otherwise
        error('Platform not defined');
end
% Object files
disp('Compiling: modnaff.c');
%mex LDFLAGS='-pthread -shared -m64' -I/usr/local/matlab/extern/include -O -c modnaff.c

eval(['mex ', PLATFORMOPTION, '-O -c modnaff.c ']);


disp('Compiling: example.c');
%mex LDFLAGS='-pthread -shared -m64' -I/usr/local/matlab/extern/include -O -c complexe.c
eval(['mex ', PLATFORMOPTION, '-O -c modnaff.c ']);
eval(['mex ', PLATFORMOPTION, '-O -c complexe.c ']);

disp('Compiling: nafflib.c');

switch computer
    case {'MACI64', 'GLNX86', 'GLNX64', 'GLNXA64'}
        internal_cc('nafflib.c modnaff.o complexe.o');
    case {'PCWIN', 'PCWIN64'}
        internal_cc('nafflib.c modnaff.obj complexe.obj');
end

cd(cd_old);

function internal_cc(fn)
% cc(filename)
%
% MAC 64 bits 
% TODO WINDOWS

disp(['Compiling: ',fn]);

switch computer
    case {'GLNX86', 'GLNX64', 'GLNXA64'}
        cmdstr = [ 'mex -I' matlabroot '/extern/include  -fPIC -O ', fn ];        
    case {'MACI64'}
        cmdstr = [ 'mex -I' matlabroot '/extern/include -O ', fn ];
    case {'PCWIN', 'PCWIN64'}
        cmdstr = [ 'mex -I' 'LDFLAGS=''-pthread -shared -m64''  -O ', fn];
    otherwise
        error('Architecture not defined')
end
disp(cmdstr);
eval(cmdstr);
