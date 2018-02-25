function [rel,tel,trand] = bpm_matrices(bpms)
%BPM_MATRICES Computes the transformation matrices for the given BPMs
%
%  INPUTS
%  1. bpms BPM indices
%
%  OUTPUTS
%  1. rel
%  2. tel
%  3. trand
%
%  EXAMPLES
%  1. [R,T,TRAND]=BPM_MATRICES(BPMS)

%
%% by L.F. Jun 2016, ESRF,  K:\machine\matlab\atlf


[rel,tel,trand]=cellfun(@extract,bpms','UniformOutput',false);

    function [rel,tel,trand]=extract(el)
        if isfield(el,'R1')
            r1=el.R1([1 3],[1 3]);
        else
            r1=eye(2);
        end
        if isfield(el,'T1')
            t1=el.T1([1 3]);
        else
            t1=[0;0];
        end
        if isfield(el,'Offset')
            tb=el.Offset(:);
        else
            tb=[0;0];
        end
        if isfield(el,'Rotation')
            C=cos(el.Rotation(1));
            S=sin(el.Rotation(1));
            rb=[C S;-S C];
        else
            rb=eye(2);
        end
        if isfield(el,'Scale')
            scale=el.Scale(:);
        else
            scale=[1;1];
        end
        if isfield(el,'Reading')
            trand=el.Reading(:);
        else
            trand=[0;0];
        end
        tel=scale.*(rb*r1*t1+tb);
        rel=[scale scale].*(rb*r1);
    end
end
