function files = passmethods(varargin);
%PASSMETHODS returns a list of available AT passmethod functions in
% /simulator/element directory
% LIST = PASSMETHODS('mex') returns list of mex-files (default)
% LIST = PASSMETHODS('m')   returns list of m-files


if nargin 
    t = varargin{1};
    if ~any(strcmpi(t,{'m','mex'}))
        t = 'mex';
    end
else
    t = 'mex';
end

switch(t)
    case 'm'
        WC = '.m';
        
    case 'mex'
        WC = ['.',mexext];
        
end
        
    
passmethoddir = fileparts(mfilename('fullpath'));
files = cellstr(ls(fullfile(passmethoddir,['*Pass',WC])));
files = strrep(files,WC,'');


