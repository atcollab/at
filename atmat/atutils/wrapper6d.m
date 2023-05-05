function varargout = wrapper6d(ring,func,varargin)
%WRAPPER6D  Private. Handle off-momentum for 6D lattice
%
%   VARARGOUT=WRAPPER6D(RING,FUNC,VARARGIN)
%
%FUNC   Wrapped function, called as FUNC(RING,IS6D,VARARGIN{:})

[warningdp6d,varargs]=getoption(varargin,'WarningDp6D');
is6d=check_6d(ring);
if is6d
    [dpargs,varargs]=getoption(varargs,{'dp','dct','df'});
    if ~isempty(dpargs)
        if warningdp6d
            warning('AT:Dp6D','\n%s\n%s\n%s',...
                'Specifying "dp" for a 6D lattice creates a copy with a modified RF frequency.',...
                'For a better efficiency, handle the RF frequency beforehand,',...
                'or to avoid this warning, use "setoption(''WarningDp6D'',false)"');
        end
        [cavargs,~]=getoption(varargs,{'cavpts'});
        ring2=atsetcavity(ring,'Frequency','nominal',dpargs{:},cavargs{:});
        [varargout{1:nargout}]=func(ring2,is6d,varargs{:});
    else
        [varargout{1:nargout}]=func(ring,is6d,varargs{:});
    end
else
    [varargout{1:nargout}]=func(ring,is6d,varargs{:});
end
end