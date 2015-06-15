function ringapert=SetPhysicalAperture(ring,apertureX,apertureY)
%ringapert=SetPhysicalAperture(ring,apertureX,apertureY)
%
% defines an aperture element before and after each element in ring
% describing the chamber dimensions.
% 
% apertureX,apertureY are the chamber vertical and horizontal aperture.
%
% chamber size will be [-apertureX,apertureX] in the horizontal plane  
% chamber size will be [-apertureY,apertureY] in the vertial plane  
%
% length(apertureX)==length(apertureY)==length(ring)
% 
% Example:
% define apertures.
% Xapert=0.06*ones(size(ring));
% Yapert=0.04*ones(size(ring));
% ring=SetPhysicalAperture(ring,Xapert/2,Yapert/2);
%
% atplot(ringinj,@plotAperture);
%
%See also: ATAPERTURE, plotAperture


if length(apertureX)==length(ring) &&  length(apertureY)==length(ring)
    ap=ataperture('AP',...
        mat2cell([-apertureX,+apertureX,-apertureY,+apertureY]...
        ,ones(size(ring)),4));
    
    ap=mat2cell(ap,ones(size(ring)));
    
    %ringapert=reshape([ap';ring';ap'],length(ring)*3,1);
    ringapert=reshape([ap';ring'],length(ring)*2,1); % only one aperture restriction
else
    error('size of ring and apertureX or apertureY are not the same')
end
    
return



