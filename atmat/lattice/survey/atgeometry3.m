function  posdata = atgeometry3(ring,varargin)
%ATGEOMETRY3 Computes the 3-D position of all elements
%
%POSDATA=ATGEOMETRY3(RING,REFPTS)
%
%RING:	 AT structure
%REFPTS: observation points (array of indexes or logical mask)
%        The allowed range is 1 to length(RING)+1
%        Defaults to 1:length(RING)+1
%
%POSDATA:Structure array, same length as REFPTS, with 3 fields:
%           x, y, z
%
%POSDATA=ATGEOMETRY3(RING,REFPTS,OFFSET)
%       Adds OFFSET(1) to the x position.
%       OFFSET(2) to the y position, OFFSET(3) to the z position
%       a scalar offset value is equivalent to [0 OFFSET 0]
%
%See also: ATGEOMETRY

[refpts,offset]=parseargs({1:length(ring)+1,[0 0 0]},varargin);
conv=eye(3);
xyzc=[0;0;0];
[xx,yy,zz]=cellfun(@incr,ring);
if isscalar(offset)
    offset=[0 offset 0];
end
xx=[0;xx]+offset(1);
yy=[0;yy]+offset(2);
zz=[0;zz]+offset(3);
posdata=struct('x',num2cell(xx(refpts)),'y',num2cell(yy(refpts)),...
    'z',num2cell(zz(refpts)));

    function varargout=incr(elem)
        if isfield(elem,'R1')
            rots(elem.R1);
        end
        if isfield(elem,'BendingAngle')
            ang=0.5*elem.BendingAngle;
            L=elem.Length/ang*sin(ang);
            hkick(-ang);
            xyzc=xyzc+conv*[L;0;0];
            hkick(-ang);
        elseif isfield(elem,'Length')
            L=elem.Length;
            xyzc=xyzc+conv*[L;0;0];
        end
        if isfield(elem,'R2')
            rots(elem.R2);
        end
        varargout=num2cell(xyzc');
    end
    function rots(R)
        cns=R(1,1);
        sns=R(1,3);
        conv=conv*[1 0 0;0 cns -sns;0 sns cns];
    end
    function hkick(ang)
        cns=cos(ang);
        sns=sin(ang);
        conv=conv*[cns -sns 0;sns cns 0;0 0 1];
    end
end

