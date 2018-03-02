function elem=ataperture(fname,varargin)
%ATAPERTURE Creates a aperture element
% ATAPERTURE(FAMNAME,LIMITS,PASSMETHOD
% define physical aperture element (collimator)
% lim=[-x,+x,-y,+y];
% 
% lim={[-x,+x,-y,+y],[-x,+x,-y,+y],[-x,+x,-y,+y],...};
% will generate various aperture elements (one for every set of errors)
% 
%  See also setphysicalaperture
% 

[rsrc,limits,method,~]=decodeatargs({[0 0],'AperturePass',''},varargin);
[limits,rsrc]=getoption(rsrc,'Limits',limits);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','Aperture');
elem=atbaselem(fname,method,'Class',cl,'Limits',limits,rsrc{:});
end
