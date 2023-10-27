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
%       Note: this is different from the radius usually defined as
%       circumference/2/pi
%
%POSDATA=ATGEOMETRY(...,'centered')
%       The offset is set so that the origin is at the centre of the ring
%
%POSDATA=ATGEOMETRY(RING,REFPTS,OFFSET)
%       Start at x=offset(1), y=offset(2). Ignored if 'centered' is set.
%       A scalar offset value is equivalent to [0 OFFSET].
%
%POSDATA=ATGEOMETRY(...,'Hangle',h_angle)
%       Set the initial trajectory angle
%
%
%See also: ATGEOMETRY3

epsil=1.e-3;
[centered,args]=getflag(varargin,'centered');
[theta0,args]=getoption(args,'Hangle',0);
[refpts,offset]=getargs(args,1:length(ring)+1,[0 0]);
xc=0;
yc=0;
thetac=theta0;
if isscalar(offset)
    offset=[0 offset];
end
[xx,yy,txy,lg,tr]=cellfun(@incr,[{struct()};ring]);    % Add a dummy element to get the origin
dff=mod(thetac+epsil, 2*pi)-epsil;
if abs(dff) < epsil         % Full ring
    xcenter=mean(xx);
    ycenter=mean(yy);
elseif abs(dff-pi) < epsil  % Half ring
    xcenter=0.5*xc;
    ycenter=0.5*yc;
else                        % Single cell
    num=cos(thetac)*xc + sin(thetac)*yc;
    den=sin(thetac-theta0);
    xcenter=-num*sin(theta0)/den;
    ycenter= num*cos(theta0)/den;
end
radius=sqrt(xcenter*xcenter+ycenter*ycenter);
if centered
    xx=xx-xcenter;
    yy=yy-ycenter;
else
    xx=xx+offset(1);
    yy=yy+offset(2);
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
