function [rel,tel,trand] = bpm_matrices(bpms)
%BPM_MATRICES Generate transformation matrices for BPM readings
%
%[REL,TEL,TRAND]=BPM_MATRICES(BPMS)
%
%Generate transformation matrices to introduce BPM errors in simulations.
%
%Positioning errors are introduced by fields of the BPM element:
%T1 and T2=-T1:  BPM positioning error
%R1 and R2=R1^1: BPM rotation
%Note: Since particle coordinates are accessed after the exit of the element,
%position errors introduced by T1,T2 are not visible since the reference orbit
%is back to nominal. This function takes care of that.
%
%Systematic BPM errors are introduced by the fields:
%Offset:         2x1 vector of H and V offsets
%Rotation        BPM rotation angle
%Scale:          2x1 vector of H and V scaling factors (calibration errors)
%
%Random errors are introduced by:
%Reading:        2x1 vector of H and V standard deviations of reading offset
%
%INPUTS:
%BPMS   Nx1 cell array of BPM elements
%
%OUTPUTS:
%REL  	1xN cell array of 2x2 rotation matrices
%TEL    1xN cell array of 2x1 translation matrices
%TRAND  1xN cell array of 2x1 random translation matrices
%
%The BPM reading vector X is obtained from the particle coordinates R by:
%
%              X = REL*R([1 3]) + TEL + TRAND.*RANDN(2,1);
%
%See also: bpm_process

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
        tel=scale.*(rb*r1*t1(:)+tb(:));
        rel=[scale scale].*(rb*r1);
    end
end
