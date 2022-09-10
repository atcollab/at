function [is_6d,ring]=check_6d(ring,enable,varargin)
%CHECK_6D	Check the presence of longitudinal motion in a lattice
%
%IS_6D=CHECK_6D(RING) Return the radiation state of RING (true/false).
%   Equivalent to IS_6D=<a href="matlab:help atGetRingProperties">atGetRingProperties</a>(RING,'is_6d')
%
%IS_6D=CHECK_6D(RING,ENABLE) Generate an error if IS_6D is different of ENABLE
%
% ENABLE:    Desired 6D state
%
%[IS_6D,NEWRING]=CHECK_6D(RING,ENABLE,'force')
%   Convert RING to the desired state using default options
%
% ENABLE:    Desired 6D state
%
% See also: ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

force=getflag(varargin,'force');
is_6d=atGetRingProperties(ring,'is_6d');

if nargin > 1
    if xor(is_6d,enable)
        if force
            if enable
                ring=atenable_6d(ring);
            else
                ring=atdisable_6d(ring);
            end
        else
            error('AT:Radiation',['''is_6d'' must be ' boolstring(enable)])
        end
    end
end

    function bstr=boolstring(enable)
        if enable
            bstr='true';
        else
            bstr='false';
        end
    end
end
