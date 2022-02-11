function varargout = wrapper6d(ring,func,varargin)
%WRAPPER6D  Private. Handle off-momentum for 6D lattice
%
%   VARARGOUT=WRAPPER6D(RING,FUNC,VARARGIN)
%
%FUNC   Wrapped function, called as FUNC(RING,IS6D,VARARGIN{:})

is6d=check_radiation(ring);
if is6d
    [cavpts,varargs]=getcavptsarg(varargin,ring);
    [dpargs,varargs]=getoption(varargs,{'dp','dct','df'});
    if ~isempty(dpargs)
        if getoption('WarningDp6D')
            warning('AT:Dp6D','\n%s\n%s\n%s',...
                'Specifying "dp" for a 6D lattice creates a copy with a modified RF frequency.',...
                'For a better efficiency, handle the RF frequency beforehand,',...
                'or to avoid this warning, use "setoption(''WarningDp6D'',false)"');
        end
        ring2=atsetcavity(ring,'Frequency','nominal',dpargs{:},'cavpts',cavpts);
        [varargout{1:nargout}]=func(ring2,is6d,varargs{:},'cavpts',cavpts);
    else
        [varargout{1:nargout}]=func(ring,is6d,varargs{:},'cavpts',cavpts);
    end
else
    [varargout{1:nargout}]=func(ring,is6d,varargin{:});
end
end