function rootdir = atroot
%ATROOT returns Accelerator Toolbox root directory
[pdir,pname,pext,pvers]=fileparts(which(mfilename));
rootdir=fullfile(pdir,'');
