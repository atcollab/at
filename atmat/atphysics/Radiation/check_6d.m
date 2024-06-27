function [is_6d, ring] = check_6d(ring,enable,varargin)
%CHECK_6D Checks the presence of longitudinal motion in a lattice.
%
%IS_6D = CHECK_6D(RING) Returns the radiation state of RING (true/false).
%  Equivalent to IS_6D=<a href="matlab:help atGetRingProperties">atGetRingProperties</a>(RING,'is_6d')
%
%IS_6D = CHECK_6D(RING, ENABLE)
%  Generates an error if IS_6D is different of ENABLE
%
%[IS_6D, NEWRING] = CHECK_6D(RING,ENABLE,'force')
%  The keyword 'force' overrides any error check, and converts
%  the RING according to the value of ENABLE.
%  IS_6D contains the status of the RING before conversion.
%
%[IS_6D, NEWRING] = CHECK_6D(RING,ENABLE,'strict',STRICT)
%  Default, STRICT=true.
%  If STRICT is true, a difference btw IS_6D and ENABLE produces an error.
%  If STRICT is false, a difference btw IS_6D and ENABLE produces a warning,
%  and the RING is converted according to the value of ENABLE.
%  IS_6D contains the status of the RING before conversion.
%
% See also: ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

[force, varargs] = getflag(varargin, 'force');
[strict, ~] = getoption(varargs, 'strict', true);
is_6d = atGetRingProperties(ring, 'is_6d');

if nargin == 1, return; end

if ~xor(is_6d,enable), return; end

if ~force
    outmsg = join(["is_6d must be" boolstring(enable) newline]);
    outerrmsg = {'AT:Radiation', outmsg};
    outwarnmsg =  {'AT:Radiation', ...
        join([outmsg ...
            "This will throw an error in a future version." newline ...
            "Set warning('off', 'AT:Radiation')" newline ...
            "or use the 'force' option to avoid this message." ...
            ]) ...
        };
    if strict, error(outerrmsg{:}); else, warning(outwarnmsg{:}); end
end

if enable, atenable_6d(ring);
else, ring = atdisable_6d(ring); end


function bstr=boolstring(enable)
    if enable
        bstr='true';
    else
        bstr='false';
    end
end
end