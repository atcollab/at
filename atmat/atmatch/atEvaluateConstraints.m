function Val=atEvaluateConstraints(R,evalfunc,posarray,indinposarray,twissin,DPP)
% This funciton evaluates the contraints defined in Constraints for lattice
% THERING
%
% Constraints: structure array struct(...
%                     'Fun',@functname(ring,lindata,globaldata),
%                     'Min',min,
%                     'Max',max,
%                     'Weight',w,
%                     'RefPoints',refpts);
%
% lindata is the output of atlinopt at the requested locations
% globdata.fractune=tune fromk atlinopt
% globdata.chromaticity=chrom from atlinopt
%
% functname must return a row vector of values to be optimized
%
% min, max and weight must have the same size as the return value of
% functname
%

% History of changes
% created 26-8-2012
% updated 28-8-2012 vector res handling.
% updated 28-8-2012 add Weigth field.
% updated 11-10-2012 any number of additonal parameters in cell array format.
%                    if only one parameter still ok.
% updated 11-03-2013 proper use of anonytmous functions introduced.
% updated 21-03-2013 Constraints now is a structure array.
%                    Weigth divides instead of multipling a constraint.
% updated 25-03-2013 Constraints call optimized for linear optics.
%

if ~isempty(posarray)
    if isempty(twissin) % recursive solution
        [lindata,t,c]=atlinopt(R,DPP,posarray);
        globdata.fractune=t;
        globdata.chromaticity=c;
    else % open line
        [lindata]=twissline(R,DPP,twissin,posarray,'chrom');
        dpp=0.0001;
        [lindatadpp]=twissline(R,DPP+dpp,twissin,posarray,'chrom');
        [lindatadpm]=twissline(R,DPP-dpp,twissin,posarray,'chrom');
        globdata.fractune=lindata(end).mu/2/pi;
        globdata.chromaticity=(lindatadpp(end).mu/2/pi-lindatadpm(end).mu/2/pi)/2/dpp;
        
    end
else
    lindata=[];
    globdata=[];
end

Val=cellfun(@(func,refpts) func(R,lindata(refpts),globdata),...
    evalfunc,indinposarray,'UniformOutput',false);
end
