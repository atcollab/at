function elem = attransformelem(elem,varargin)
%ATTRANSFORMELEM Set the translations and rotations of an element
%
%NEWELEM=ATTRANSFORMELEM(ELEM,KEY,VALUE,...)
%   Set the translations and rotations of ELEM, and compute the R1, T1, R2
%   and T2 fields from them.
%
%KEYS:
%   'dx'            Horizontal shift [m]. Default: no change
%   'dy'            Vertical shift [m]. Default: no change
%   'dz'            Longitudinal shift [m]. Default: no change
%   'tilt'          Tilt angle, rotation around the z-axis [rad]. Default: no change
%   'pitch'         Pitch angle, rotation around the x-axis [rad]. Default: no change
%   'yaw'           Yaw angle, rotation around the y-axis [rad]. Default: no change
%   'tilt_frame'    Tilt angle of the reference frame [rad]. Useful to
%                   generate vertical bending magnets. Default: no change
%   'reference'     Transformation reference, 'centre' or 'entrance'.
%                   Default: the reference already stored in the element,
%                   otherwise 'centre'
%
%NEWELEM=ATTRANSFORMELEM(...,'relative')
%   The input values are added to the previous ones.
%
%NEWELEM=ATTRANSFORMELEM(ELEM)
%   Recompute R1, T1, R2 and T2 from the transformation fields stored in
%   ELEM. Use this after editing these fields directly.
%
%   The axes follow the right hand rule: x is horizontal, y is vertical and
%   z is along the beam. A positive angle is a clockwise rotation when
%   looking in the direction of the rotation axis.
%
%   The translations are applied before the rotations. The rotations are
%   applied in the order tilt -> yaw -> pitch. The element is rotated around
%   its entrance or its centre (the middle of the chord joining the entry and
%   exit points of the element).
%
%   The transformation is stored in the fields dx, dy, dz, tilt, pitch, yaw,
%   tilt_frame and reference (0 for centre, 1 for entrance), with the same
%   names and meaning as in pyAT. R1, T1, R2 and T2 are recomputed from these
%   fields, any previous value is discarded.
%
%   The transverse momenta are canonical, so the angular kicks of pitch and
%   yaw scale with (1 + delta). This introduces spurious dispersion.
%
%   The implementation follows doi:10.1016/j.nima.2022.167487 and pyAT's
%   transform_elem. Comments featuring 'Eq' point to the paper's equations.
%
%See also atsettransform, atshiftelem, attiltelem

names={'dx','dy','dz','pitch','yaw','tilt','tilt_frame'};
[relative,args]=getflag(varargin,'relative');
[reference,args]=getoption(args,'reference',[]);
newvals=cell(size(names));
for i=1:length(names)
    [newvals{i},args]=getoption(args,names{i},[]);
end
if ~isempty(args)
    error('AT:WrongArgument','Unexpected argument "%s"',num2str(args{1}))
end

if isfield(elem,'reference')
    reference0=elem.reference;
else
    reference0=[];
end
if isempty(reference)
    if isempty(reference0)
        reference=0;
    else
        reference=reference0;
    end
else
    reference=refcode(reference);
end
if relative && ~isempty(reference0) && reference ~= reference0
    error('AT:WrongArgument',['Element %s: reference point change ',...
        'not allowed for relative transformations'],elem.FamName);
end

vals=zeros(size(names));
for i=1:length(names)
    if isfield(elem,names{i})
        vals(i)=elem.(names{i});
    end
    if ~isempty(newvals{i})
        if relative
            vals(i)=vals(i)+newvals{i};
        else
            vals(i)=newvals{i};
        end
    end
    elem.(names{i})=vals(i);
end
elem.reference=reference;
offsets=vals(1:3)';
rotations=vals(4:6);
tilt_frame=vals(7);

zaxis=[0;0;1];
if isfield(elem,'Length'), L=elem.Length; else, L=0; end
if isfield(elem,'BendingAngle'), angle=elem.BendingAngle; else, angle=0; end

% Rotated transfer matrix linked to a bending angle
RB=rotation(0,-angle,0);            % Eq. (12)
RB_half=rotation(0,-angle/2,0);     % Eq. (12)

if reference == 0
    % Compute entrance rotation matrix in the rotated frame
    r3d_entrance=RB_half*rotation(rotations(1),rotations(2),rotations(3))*RB_half';  % Eq. (31)
    if angle ~= 0
        Rc=L/angle;
        OO0=Rc*sin(angle/2)*RB_half*zaxis;                      % Eq. (34)
        P0P=-Rc*sin(angle/2)*r3d_entrance*RB_half*zaxis;        % Eq. (36)
    else
        OO0=L/2*zaxis;                                          % Eq. (34)
        P0P=-L/2*r3d_entrance*zaxis;                            % Eq. (36)
    end
    % Transform offset to magnet entrance
    OP=OO0+P0P+RB_half*offsets;                                 % Eq. (33)
else
    r3d_entrance=rotation(rotations(1),rotations(2),rotations(3));   % Eq. (3)
    OP=offsets;                                                 % Eq. (2)
end

% R1, T1
ld_entrance=r3d_entrance(:,3)'*OP;                              % Eq. (33)
R1=r_matrix(ld_entrance,r3d_entrance);
T1=R1\translation_vector(ld_entrance,r3d_entrance,OP,r3d_entrance(:,1),r3d_entrance(:,2));

% R2, T2
r3d_exit=RB'*r3d_entrance'*RB;                                  % Eq. (18) or (32)
if angle ~= 0
    Rc=L/angle;
    OPp=[Rc*(cos(angle)-1); 0; L*sin(angle)/angle];             % Eq. (24)
else
    OPp=L*zaxis;                                                % Eq. (24)
end
OOp=r3d_entrance*OPp+OP;                                        % Eq. (25)
OpPp=OPp-OOp;
ld_exit=RB(:,3)'*OpPp;                                          % Eq. (23) or (37)
R2=r_matrix(ld_exit,r3d_exit);
T2=translation_vector(ld_exit,r3d_exit,OpPp,RB(:,1),RB(:,2));

tf=tilt_frame_matrix(tilt_frame);
elem.R1=tf*R1;
elem.R2=R2*tf';
elem.T1=T1;
elem.T2=T2;

    function code=refcode(ref)
        if ischar(ref) || isstring(ref)
            switch lower(char(ref))
                case {'centre','center'}
                    code=0;
                case 'entrance'
                    code=1;
                otherwise
                    error('AT:WrongArgument',['Unsupported reference "%s", ',...
                        'please choose either ''centre'' or ''entrance''.'],ref);
            end
        elseif isscalar(ref) && (ref == 0 || ref == 1)
            code=double(ref);
        else
            error('AT:WrongArgument',...
                'Unsupported reference, please choose either ''centre'' or ''entrance''.');
        end
    end
end

function r=rotation(alpha,beta,gamma)
% 3D rotation matrix using the Tait-Bryan angles convention, Eq. (3)
% alpha: rotation about the x-axis (pitch)
% beta: rotation about the y-axis (yaw)
% gamma: rotation about the z-axis (tilt)
Rx=[1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Ry=[cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz=[cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
r=Rx*Ry*Rz;
end

function t=translation_vector(ld,r3d,offsets,xaxis,yaxis)
% Translation vector resulting from the joint effect of a longitudinal
% displacement (in the rotated frame), 3D offsets, and the 3D rotation
% matrix, Eqs. (8-11)
% ld: longitudinal displacement [m]
% xaxis, yaxis: unit axes of the rotated frame in the xyz coordinate system
tD0=[-offsets'*xaxis; 0; -offsets'*yaxis; 0; 0; 0];
T0=[ld*r3d(3,1)/r3d(3,3); r3d(3,1); ld*r3d(3,2)/r3d(3,3); r3d(3,2); 0; ld/r3d(3,3)];
t=T0+tD0;
end

function m=r_matrix(ld,r3d)
% Rotation matrix operator (R1, R2), including the effect of a longitudinal
% displacement ld, Eq. (9)
c2=r3d(3,3)^2;
m=[r3d(2,2)/r3d(3,3) ld*r3d(2,2)/c2 -r3d(1,2)/r3d(3,3) -ld*r3d(1,2)/c2 0 0;...
    0 r3d(1,1) 0 r3d(2,1) r3d(3,1) 0;...
    -r3d(2,1)/r3d(3,3) -ld*r3d(2,1)/c2 r3d(1,1)/r3d(3,3) ld*r3d(1,1)/c2 0 0;...
    0 r3d(1,2) 0 r3d(2,2) r3d(3,2) 0;...
    0 0 0 0 1 0;...
    -r3d(1,3)/r3d(3,3) -ld*r3d(1,3)/c2 -r3d(2,3)/r3d(3,3) -ld*r3d(2,3)/c2 0 1];
end

function rm=tilt_frame_matrix(rots)
cs=cos(rots);
sn=sin(rots);
rm=diag([cs cs cs cs 1 1]);
rm(1,3)=sn;
rm(2,4)=sn;
rm(3,1)=-sn;
rm(4,2)=-sn;
end
