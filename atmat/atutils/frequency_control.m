function varargout = frequency_control(func,ring,varargin)
%FREQUENCY_CONTROL  Private. Handle off-momentum for 6D lattice
%
%   VARARGOUT=FREQUENCY_CONTROL(FUNC,RING,VARARGIN)
%
%FUNC   Wrapped function, called as FUNC(RING,VARARGIN{:},'is6d',IS6D)

[warningdp6d,varargs]=getoption(varargin,'WarningDp6D');
[is_6d, varargs]= getoption(varargs, 'is_6d', []);
if isempty(is_6d), is_6d=atGetRingProperties(ring,'is_6d'); end
if is_6d
    [dpargs,varargs]=getoption(varargs,{'dp','dct','df'});
    if length(dpargs) > 2
        error('AT:OffMomentum', ['Off-momentum specification: ',...
            'only one keyword in "dp", "dct" or "df" may be specified.'])
    end
    if ~isempty(dpargs)
        if warningdp6d
            warning('AT:Dp6D','\n%s\n%s\n%s',...
                'Specifying "dp, dct or df" for a 6D lattice creates a copy with a modified RF frequency.',...
                'For a better efficiency, handle the RF frequency beforehand,',...
                'or to avoid this warning, use "setoption(''WarningDp6D'',false)"');
        end
        [cavargs,~]=getoption(varargs,{'cavpts'});
        ring2=atsetcavity(ring,'Frequency','nominal',dpargs{:},cavargs{:});
        [varargout{1:nargout}]=func(ring2,varargs{:},'is_6d',is_6d);
    else
        [varargout{1:nargout}]=func(ring,varargs{:},'is_6d',is_6d);
    end
else
    [varargout{1:nargout}]=func(ring,varargs{:},'is_6d',is_6d);
end
end