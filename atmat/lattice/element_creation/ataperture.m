function elem=ataperture(fname,varargin)
% define physical aperture element (collimator)
% lim=[-x,+x,-y,+y];
% 
% lim={[-x,+x,-y,+y],[-x,+x,-y,+y],[-x,+x,-y,+y],...};
% will generate various aperture elements (one for every set of errors)
% 
%See also: SETPHYSICALAPERTURE
% 

[rsrc,limits,method,~]=decodeatargs({[0 0],'AperturePass',''},varargin);
[rsrc,limits]=getatarg(rsrc,limits,'Limits');
elem=atbaselem(fname,method,'Limits',limits,'Class','Aperture',rsrc{:});
end