function rootdir = atroot
%ATROOT returns Accelerator Toolbox root directory

w = which('DriftPass');
if ~isempty(w)
    d = fileparts(w);
    fileparts(w);
    ind = strfind(d,[filesep,'simulator', filesep, 'element'])-1;
    rootdir = d(1:ind);
else
    error('Unable to find AT root directory - DriftPass.m not on MATLAB path');
end

    
