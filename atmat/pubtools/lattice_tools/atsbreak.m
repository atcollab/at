function [newring,refpts] = atsbreak(ring,sdata)
%ATSBREAK Insert markers at given s positions in a lattice
%
%[NEWRING,REFPTS]=ATSBREAK(RING,SPOS) build a new lattice with breakpoints
%
%RING:  AT structure
%SPOS:  Poition of the desired breakpoints
%
%NEWRING:   Modified AT structure
%REFPTS:    Index of breakpoints
%
% See also ATINSERTELEMS ATSPLITELEM ATDIVELEM

mk=atmarker('xyzt');
selem=findspos(ring,1:length(ring)+1);
[sv,isort]=sort(mod(sdata(:)',selem(end)));
ir=1;
is=1;
it=0;
nr={};
for ie=1:length(ring)
    if selem(ie+1) > sv(is)
        elem=ring{ie};
        for j=is:length(sv)
            if selem(ie+1) <= sv(j)
                last=j-1;
                break
            end
            last=j;
        end
        frac=(sv(is:last)-selem(ie))./elem.Length;
        is=last+1;
        it=it+1;
        nr{it}=[ring(ir:ie-1);atsplitelem(elem,frac,mk)]; %#ok<AGROW>
        ir=ie+1;
        if is > length(sv), break; end
    end
end
newring=cat(1,nr{:},ring(ir:end));
refpts(isort)=find(atgetcells(newring,'FamName','xyzt'));
refpts=reshape(refpts,size(sdata));
end
