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
        ring2=atsetcavity(ring,'Frequency','nominal',dpargs{:},'cavpts',cavpts);
        [varargout{1:nargout}]=func(ring2,is6d,varargs{:},'cavpts',cavpts);
    else
        [varargout{1:nargout}]=func(ring,is6d,varargs{:},'cavpts',cavpts);
    end
else
    [varargout{1:nargout}]=func(ring,is6d,varargin{:});
end
end