function [value,opts] = getoption(opts,key,value)
%GETOPTION Extract a keyword argument from an argument list
%
%VALUE=GETOPTION(ARGS,'KEY',DEFAULT)
%VALUE=GETOPTION(ARGS,KEY=DEFAULT)  in Matlab >= R2021a
%   Extract a keyword argument, in the form of a pair "key,value" from
%   input arguments ARGS (typically from VARARGIN).
%   Return DEFAULT value if the keyword is absent
%
% ARGS:     Argument list: cell array (usually VARARGIN) or structure
% KEY:      Key name
% DEFAULT:  Value used if "key,value" is absent from the argument list
%
%VALUE=GETOPTION(ARGS,'KEY')
%   The default value is taken from a list of predefined keys. Use
%   GETOPTION() for the list of predefined keys
%
%VALUE=GETOPTION(ARGS,{'KEY1','KEY2',...)
%   Value is the list of key/value pairs matching KEY1 or KEY2 or...
%
%VALUE=GETOPTION('KEY')
%   Return the default value of a predefined key. Use GETOPTION() for
%   the list of predefined keys
%
%VALUE=GETOPTION()
%   Return all the default values
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
%See also GETFLAG, GETARGS, SETOPTION, ATOPTIONS

if nargin < 1           % list of all available predefined keys
    value=atoptions.instance();
elseif nargin < 2       % value of the predefined key
    value=atoptions.instance().(opts);
elseif iscell(key)      % extract a subset of key/value pairs
    oksel=false(size(opts));
    for k=key
        ok=[cellfun(@(v) (ischar(v) || isstring(v)) && strcmpi(v, k{1}), ...
            opts(1:end-1)) false];  % option name cannot be the last argument
        oksel=oksel | ok | circshift(ok,[0,1]);
    end
    value=opts(oksel);
    opts(oksel)=[];
else                    % return the value for the given key
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
