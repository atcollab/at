function slicedring = atslice(ring,npts)
%ATSLICE Splits ring elements in short slices for plotting

totlength=findspos(ring,length(ring)+1);
elmlength=totlength/npts;
slicedring={};
for i=1:length(ring)
    elem=ring{i};
    if isfield(elem,'BendingAngle')
        newelems=splitdipole(elem,elmlength);
    elseif elem.Length > 0
        nslices=ceil(elem.Length/elmlength);
        elem.Length=elem.Length/nslices;
        elt={elem};
        newelems=elt(ones(nslices,1));
    else
        newelems={elem};
    end
    slicedring=[slicedring;newelems]; %#ok<AGROW>
end

%   Special treatment of dipoles

    function newelems=splitdipole(elem,elmlength)
        nsl=ceil(elem.Length/elmlength);
        if isfield(elem,'EntranceAngle')
            ena=elem.EntranceAngle;
            elem.EntranceAngle=0;
        end
        if isfield(elem,'ExitAngle')
            exa=elem.ExitAngle;
            elem.ExitAngle=0;
        end
        elem.Length=elem.Length/nsl;
        elem.BendingAngle=elem.BendingAngle/nsl;
        el={elem};
        newelems=el(ones(nsl,1));
        if exist('ena','var'), newelems{1}.EntranceAngle=ena; end
        if exist('exa','var'), newelems{end}.ExitAngle=exa; end
    end
end

