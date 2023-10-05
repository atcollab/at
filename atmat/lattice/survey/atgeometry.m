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
%       Start at x=offset(1), y=offset(2)
%       a scalar offset value is equivalent to [0 OFFSET]
%
%POSDATA=ATGEOMETRY(...,'centered')
%       The offset is set as [0 RADIUS]
%
%POSDATA=ATGEOMETRY(...,'Hangle',h_angle)
%       Set the initial trajectory angle
%
%
%See also: ATGEOMETRY3

[centered,args]=getflag(varargin,'centered');
[theta0,args]=getoption(args,'Hangle',0);
[refpts,offset]=getargs(args,1:length(ring)+1,[0 0]);
xc=0;
yc=0;
thetac=theta0;
if isscalar(offset)
    yc=offset;
else
    xc=offset(1);
    yc=offset(2);
end
[xx,yy,txy,lg,tr]=cellfun(@incr,[{struct()};ring]);    % Add a dummy element to get the origin
radius=-(yc+xc/tan(thetac-theta0));
if centered
    xx=xx-radius*sin(theta0);
    yy=yy+radius*cos(theta0);
end
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
        c=cos(thetac); s=sin(thetac);
        lgtr=[xc yc]*[c -s;s c];
        varargout={xc,yc,thetac,lgtr(1),lgtr(2)};
    end
end
