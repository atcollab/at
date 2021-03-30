function [radon,ring]=check_radiation(ring,onoff,varargin)
%CHECK_RADIATION	Check the radiation state of a ring
%
%RAD=CHECK_RADIATION(RING) Return the radiation state of RING (true/false)
%
%RAD=CHECK_RADIATION(RING,ONOFF) Generate an error if RAD is different of ONOFF
%
% ONOFF:    Desired radiation state
%
%[RAD,NEWRING]=CHECK_RADIATION(RING,ONOFF,'force')
%   Convert RING to the desired state using default options
%
% ONOFF:    Desired radiation state

force=getflag(varargin,'force');
radon=false;
for i=1:length(ring)
    passmethod=ring{i}.PassMethod;
    if endsWith(passmethod, {'RadPass', 'CavityPass', 'QuantDiffPass'})
        radon=true;
        break;
    end
end

if nargin >= 2
    if xor(radon,onoff)
        if force
            if onoff
                ring=atradon(ring);
            else
                ring=atradoff(ring);
            end
        else
            error('AT:Radiation',['Radiation must be ' boolstring(onoff)])
        end
    end
end

    function bstr=boolstring(onoff)
        if onoff
            bstr='ON';
        else
            bstr='OFF';
        end
    end
end
