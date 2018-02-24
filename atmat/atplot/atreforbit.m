function  [xref,zref] = atreforbit(ring)
%ATREFORBIT - computes the coordinates of the local referential
%           It allows plotting functions (trajectory/orbit) to plot through
%           displaced elements
%  INPUTS
%    1. ring Ring structure
%
%  OUTPUTS
%    1. xref Horizontal reference orbit
%    2. zref Vertical reference orbit
%
%  EXAMPLE
%    1. [XREF,ZREF]=ATREFORBIT(RING)
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
        slope=slope-tan(T1([2 4]));
    end
end

