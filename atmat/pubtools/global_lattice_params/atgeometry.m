function  [posdata,radius] = atgeometry(ring,varargin)
%ATGEOMETRY Computes the 2-D position of all elements (no vertical bend)
%
%POSDATA=ATGEOMETRY(RING,REFPTS)
%
%RING:	 AT structure
%REFPTS: observation points (array of indexes or logical mask)
%        The allowed range is 1 to length(RING)+1
%        Defaults to 1:length(RING)+1
%
%POSDATA:Structure array, same length as REFPTS, with 5 fields:
%           x, y, angle, long, trans
%
%[POSDATA,RADIUS]=ATGEOMETRY(RING,REFPTS)
%       Outputs the machine radius at the beginning of the lattice.
%       Meaningful if RING is a cell of a periodic lattice.
%
%POSDATA=ATGEOMETRY(RING,REFPTS,OFFSET)
%       Adds OFFSET(1) to the x position and OFFSET(2) to the y position
%       a scalar offset value is equivalent to [0 OFFSET]
%POSDATA=ATGEOMETRY(RING,REFPTS,'centered')
%       The offset is set as [0 RADIUS]
%
%
%See also: ATGEOMETRY3

[refpts,offset]=parseargs({1:length(ring)+1,[0 0]},varargin);
xc=0;
yc=0;
thetac=0;
[xx,yy,txy]=cellfun(@incr,ring);
radius=-(yc+xc/tan(thetac));
if ischar(offset) && strcmp(offset,'centered')
    offset=[0 radius];
elseif isscalar(offset)
    offset=[0 offset];
end
xx=[0;xx]+offset(1);
yy=[0;yy]+offset(2);
txy=[0;txy];
[lg,tr]=arrayfun(@srcpt,xx,yy,txy);
posdata=struct('x',num2cell(xx(refpts)),'y',num2cell(yy(refpts)),...
    'angle',num2cell(txy(refpts)),'long',num2cell(lg(refpts)),...
    'trans',num2cell(tr(refpts)));

    function varargout=incr(elem)
        if isfield(elem,'BendingAngle') && elem.BendingAngle ~= 0
            ang=0.5*elem.BendingAngle;
            thetac=thetac-ang;
            L=elem.Length/ang*sin(ang);
            xc=xc+L*cos(thetac);
            yc=yc+L*sin(thetac);
            thetac=thetac-ang;
        elseif isfield(elem,'Length')
            L=elem.Length;
            xc=xc+L*cos(thetac);
            yc=yc+L*sin(thetac);
        end
        varargout={xc,yc,thetac};
    end

    function varargout=srcpt(x,y,angle)
        c=cos(angle); s=sin(angle);
        varargout=num2cell([x y]*[c -s;s c]);
    end
end
