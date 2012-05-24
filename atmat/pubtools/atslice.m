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
        ena=elem.EntranceAngle;
        exa=elem.ExitAngle;
        elem.Length=elem.Length/nsl;
        elem.BendingAngle=elem.BendingAngle/nsl;
        elem.EntranceAngle=0;
        elem.ExitAngle=0;
        el={elem};
        newelems=el(ones(nsl,1));
        newelems{1}.EntranceAngle=ena;
        newelems{end}.ExitAngle=exa;
    end
end

