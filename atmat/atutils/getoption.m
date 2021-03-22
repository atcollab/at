function [value,opts] = getoption(opts,key,value)
%GETOPTION Extract a keyword argument from an argument list
%
%VALUE=GETOPTION(ARGS,'KEY',DEFAULT)
%   Extract a keyword argument, in the form of a pair "key,value" from
%   input arguments ARGS (typically from VARARGIN).
%   Return DEFAULT value if the keyword is absent
% ARGS:     Argument list: cell array (usually VARARGIN) or structure
% KEY:      Key name
% DEFAULT:  Value used if "key,value" is absent from the argument list
%
%VALUE=GETOPTION(ARGS,'KEY')
%   The default value is taken from a list of predefined keys. See
%   INITOPTIONS for the list of predefined keys
%
%VALUE=GETOPTION('KEY')
%   Return the default value of a predefined key. See INITOPTIONS for
%   the list of predefined keys
%
%[VALUE,REMARGS]=GETOPTION(ARGS,...)
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
%See also GETFLAG, GETARGS, SETOPTION, INITOPTIONS

if nargin < 1
    value=atoptions.instance();
elseif nargin < 2
    value=atoptions.instance().(opts);
else
    if nargin < 3
        value=atoptions.instance().(key);
    end
    if iscell(opts)
        ok=[cellfun(@(v) (ischar(v) || isstring(v)) && strcmpi(v, key), ...
            opts(1:end-1)) false];  % option name cannot be the last argument
        if any(ok)
            okval=circshift(ok,[0,1]);
            value=opts{find(okval,1,'last')};
            opts(ok|okval)=[];
        end
    elseif isstruct(opts)
        if isfield(opts,key)
            value=opts.(key);
            opts=rmfield(opts,key);
        end
    end
end
end
