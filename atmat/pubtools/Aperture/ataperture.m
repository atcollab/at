function ap=ataperture(fam,lim)
% define physical aperture element (collimator)
% lim=[-x,+x,-y,+y];
% 
% lim={[-x,+x,-y,+y],[-x,+x,-y,+y],[-x,+x,-y,+y],...};
% will generate various aperture elements (one for every set of errors)
% 
%See also: SETPHYSICALAPERTURE
% 

ap=struct('FamName',fam,...
    'Class','Aperture',...
    'PassMethod','AperturePass',...
    'Limits',lim,'Length',0);
