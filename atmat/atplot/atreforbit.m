function  [xref,zref] = atreforbit(ring)
%ATREFORBIT Keep track of the nominal reference orbit through displaced elements
%
%Element displacement vectors T1 and T2 are often used to introduce positioning
%errors of elements. When plotting the resulting closed orbit or trajectories,
%it is useful to plot positions with respect to the theoretical reference rather
%than to the displaced reference. This function generates reference coordinates
%xref and zref so that:
%
%- R(1) and R(3) are the particles horizontal and vertical positions with
%respect to the actual (displaced) reference,
%
%- R(1)+xref and R(3)+zref are the positions with respect to the ideal reference.
%
%If any of T1 or T2 is part of the design of the ring, the reference
%translation can be skipped by setting a field 'hideT1' or 'hideT2' on the
%element.
%
%  INPUTS
%  1. ring Ring structure
%
%  OUTPUTS
%  1. xref Horizontal reference orbit shift
%  2. zref Vertical reference orbit shift
%
%  EXAMPLE
%  1. [XREF,ZREF]=ATREFORBIT(RING)
%
%  See also atplot

xzc=[0;0];
slope=[0;0];
[xx,zz]=cellfun(@incr,ring);
xref=[0;xx];
zref=[0;zz];

    function varargout=incr(elem)
%         if isfield(elem,'R1')
%             rots(elem.R1);
%         end
        if isfield(elem,'T1') && ~isfield(elem,'hideT1')
            hvkick(elem.T1(:));
        end
        if isfield(elem,'Length')
            xzc=xzc+slope*elem.Length;
        end
        if isfield(elem,'T2') && ~isfield(elem,'hideT2')
            hvkick(elem.T2(:));
        end
%         if isfield(elem,'R2')
%             rots(elem.R2);
%         end
        varargout=num2cell(xzc');
    end
%     function rots(R)
%         cns=R(1,1);
%         sns=R(1,3);
%         conv=conv*[1 0 0;0 cns -sns;0 sns cns];
%     end
    function hvkick(T1)
        xzc=xzc-T1([1 3]);
        slope=slope-T1([2 4]);
    end
end

