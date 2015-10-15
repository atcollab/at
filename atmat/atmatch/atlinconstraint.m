function constraint = atlinconstraint(refpts,params,vmin,vmax,weight)
%ATLINCONSTRAINT Generate a constraint on linear optics for atmatch
%
%CONSTRAINT=ATLINCONSTRAINT(REFPTS,PARAMS,VMIN,VMAX,WEIGHT)
%
%REFPTS Row vector of selected positions
%PARAMS Cell array describing the desired value at each position
%       The length of params must be 1 or length(REFPTS)
%       Each element of PARAMS is itself a cell array used to select a
%       field in the structure returned by atlinopt
%VMIN   Minimum value for the constraint
%VMAX   Maximum value for the constraint
%
%CONSTRAINT Row structure array to be used in atmatch
%
% REFPTS, PARAMS, VMIN, VMAX, WEIGHT must have the same length,
%       or have length 1

ncst=max(length(refpts),length(params));
refp=expand(refpts,ncst);
pars=expand(params,ncst);
vmin=expand(vmin,ncst);
vmax=expand(vmax,ncst);
weight=expand(weight,ncst);

paramnames=cellfun(@(c) c{1},pars,'UniformOutput',false);
tunes=strcmp('tune',paramnames);
chrom=strcmp('chromaticity',paramnames);
other=~(chrom|tunes);

constraint=struct([]);
if any(tunes)
    constraint=[constraint,struct(...
        'Fun',@(~,~,gd) cellfun(@(p) gd.fractune(p{2}{:}),pars(tunes)),...
        'Min',vmin(tunes),...
        'Max',vmax(tunes),...
        'RefPoints',refp(tunes),...
        'Weight',weight(tunes))];
end
if any(chrom)
    constraint=[constraint,struct(...
        'Fun',@(~,~,gd) cellfun(@(p) gd.chromaticity(p{2}{:}),pars(chrom)),...
        'Min',vmin(chrom),...
        'Max',vmax(chrom),...
        'RefPoints',refp(chrom),...
        'Weight',weight(chrom))];
end
if any(other)
    constraint=[constraint,struct(...
        'Fun',@(~,ld,~) arrayfun(@(lindata,p) getfield(lindata,p{1}{:}),ld,pars(other)),...
        'Min',vmin(other),...
        'Max',vmax(other),...
        'RefPoints',refp(other),...
        'Weight',weight(other))];
end

    function val=expand(val,sz)
        if length(val) == 1
            val=val(ones(1,sz));
        end
    end

end
