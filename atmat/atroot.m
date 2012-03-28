function rootdir = atroot
%ATROOT returns Accelerator Toolbox root directory
[pdir,pname,pext]=fileparts(which(mfilename));
rootdir=fullfile(pdir,'');
