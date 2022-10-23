function varargout = getdparg(args,varargin)
%GETDPARG Handle positional dp arguments
%
%[DP,VARARGS]=GETDPARG(VARARGS)
%   If the 1st argument in VARARGS is a scalar numeric less than 1, it is
%   considered as DP and removed from VARARGS. 
%
%VARARGS=GETDPARG(VARARGS)
%   DP is extracted, and if it is finite and non-zero,
%   {'DP', DP} is added to VARARGS

[check,varargs]=getoption(varargin,'check',@isdparg);
[key,varargs]=getoption(varargs,'key','dp');
[value,~]=getargs(varargs,NaN);     % Default value
if ~isempty(args) && check(args{1}) % Positional dp argument
    value=args{1};
    args=args(2:end);
end
[value,args]=getoption(args,key,value);  % Keyword argument takes precedence
if nargout == 1
    if isfinite(value)
        args=[args {key, value}];
    end
    varargout={args};
else
    varargout={value,args};
end

    function ok=isdparg(arg)
        % ~(arg >= 1) : NaN is accepted as a valid dp value
        ok=isnumeric(arg) && isscalar(arg) && ~(arg >= 1);
    end
end