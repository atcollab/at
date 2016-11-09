function rerr=setGirderError(r,pert,errval,mag_group)
%        rerr=setGirderError(r,pert,errval,mag_group)
% 
%  sets specified girder errors (absolute values)
%  mag_group is the output of getMagGroupsFromGirderIndex(ring);
%
%
%see also: 
if nargin<4
    mag_group=getMagGroupsFromGirderIndex(r);
end
magindex=[mag_group{:}];

% get existing errors
[X0,Y0,~,T0,R0,P0]=GetExistingErrors(r,magindex);

% initialize T1 T2 R1 R2 errors to zero
r=setANYshift(r,magindex,1:6,zeros(size(magindex)));
r=setTiltAbout(r,magindex,'s',zeros(size(magindex)));
% initialize rotation storage field
r=atsetfieldvalues(r,magindex,'RotAboutX',{1,1},0);
r=atsetfieldvalues(r,magindex,'RotAboutY',{1,1},0);
r=atsetfieldvalues(r,magindex,'Tilt',{1,1},0);

switch pert
    case 'gx'
        errfun=@(r,po,er)setANYshift(r,po,1,er);
        err0=X0;
    case 'gy'
        errfun=@(r,po,er)setANYshift(r,po,3,er);
        err0=Y0;
    case 'gphi'
        errfun=@(r,po,er)setTiltGirderAbout(r,mag_group,'x',er);
        err0=R0;
    case 'gtheta'
        errfun=@(r,po,er)setTiltGirderAbout(r,mag_group,'y',er);
        err0=P0;
    case 'gpsi'
        errfun=@(r,po,er)setTiltGirderAbout(r,mag_group,'s',er);
        err0=T0;
end

err_group=cellfun(@(indg,ev)ev*ones(size(indg)),mag_group,num2cell(errval)','un',0);
   
errval=[err_group{:}];

% move girder
rerr=errfun(r,magindex,err0+errval);

return