function varargout = getdparg(args,varargin)
%GETDPARG Handle positional dp arguments
%
%[DP,VARARGS]=GETDPARG(VARARGS)
%   If the 1st argument in VARARGS is empty or scalar numeric, it is
%   considered as DP and removed from VARARGS. 
%
%VARARGS=GETDPARG(VARARGS)
%   DP is extracted, and if it is finite and non-zero,
%   {'DP', DP} is added to VARARGS

[check,varargs]=getoption(varargin,'check',@isdparg);
[dp,varargs]=getargs(varargs,NaN); %#ok<ASGLU> 
[dp,args]=getoption(args,'dp',dp);
[dp,args]=getargs(args,dp,'check',check);
if isempty(dp), dp=NaN; end
if nargout == 1
    if isfinite(dp)
        args=[args {'dp', dp}];
    end
    varargout={args};
else
    varargout={dp,args};
end

    function v=isdparg(arg)
        v=isempty(arg) || (isnumeric(arg) && isscalar(arg));
    end
end