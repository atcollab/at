function [value,opts] = getoption(opts,key,default)
%GETOPTION Extract a keyword argument from the input arguments
%
%VALUE=GETOPTION(ARGS,'KEY',DEFAULT)
%   Extract a keyword argument, in the form of a pair "key,value" from
%   input arguments ARGS (typically from VARARGIN).
%	Return DEFAULT value if the keyword is absent
% ARGS:     Argument list (cell array or structure, usually VARARGIN)
% KEY:      Key name
% DEFAULT:  Value when "key,value" is absent from the input arguments
%
%VALUE=GETOPTION(ARGS,'KEY')
%   The default value is taken from a list of predefined keys. See
%   ATINIPREFS for the list of predefined keys
%
%DEFAULT=GETOPTION('KEY')
%   Return the default value for a predefined key. See ATINIPREFS for
%   the list of predefined keys
%
%[VALUE,NEWARGS]=GETOPTION(ARGS,...)
%  Return the remaining arguments after removing the processed ones
%
%Example:
%
%function testfunc(varargin)
%
%[flag,args] = getflag(varargin, 'Flag');       % Extract an optional flag
%[range,args] = getoption(args, 'Range', 1:10); % Extract a keyword argument
%[width, height] = getargs(args, 210, 297});    % Extract positional arguments
%
%See also GETFLAG, GETARGS, SETOPTION, ATINIPREFS

if nargin < 2
    value=atprefutil('get',opts);
else
    if nargin < 3
        default = atprefutil('get',key);
    end
    if iscell(opts)
        ok=[strcmpi(key,opts(1:end-1)) false];  % option name cannot be the last argument
        if any(ok)
            okval=circshift(ok,[0,1]);
            default=opts{find(okval,1,'last')};
            opts(ok|okval)=[];
        end
    elseif isstruct(opts)
        if isfield(opts,key)
            default=opts.(key);
            opts=rmfield(opts,key);
        end
    end
    value=default;
end
end
