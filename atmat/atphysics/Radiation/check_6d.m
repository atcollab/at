function [is_6d, newring] = check_6d(ring,enable,varargin)
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
%  If STRICT is true, a non-6d ring produces an error.
%  If STRICT is false, a non-6d ring produces a warning,
%  and the RING is converted according to the value of ENABLE.
%  IS_6D contains the status of the RING before conversion.
%
% See also: ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

[force, varargs] = getflag(varargin, 'force');
[strict, ~] = getoption(varargs, 'strict', true);
is_6d = atGetRingProperties(ring, 'is_6d');

if nargin == 1, return; end

if xor(is_6d, enable) && strict && ~force
    error('AT:Radiation',['''is_6d'' must be ' boolstring(enable)]);
end

if ~is_6d && ~force
    if strict, error('RING is not 6d.');
    else, warning('RING is not 6d.'); end
end

if enable, newring = atenable_6d(ring);
else, newring = atdisable_6d(ring); end

function bstr=boolstring(enable)
    if enable
        bstr='true';
    else
        bstr='false';
    end
end
end