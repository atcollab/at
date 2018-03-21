function maggroups=getMagGroupsFromMagNum(r)
%GETMAGGROUPSFROMMAGNUM - Gets magnet from a Magnet group
% output maggroups in r given that the variable MagNum is available in the
% lattice.
% maggroups is a cell array of magnet indexes describing a single magnet in
% reality, but sliced in the lattice
% a single magnet has the same MagNum value.
% 
%see also SectorDipoleWithParallelFaces 

maggroups={};

magnumind=findcells(r,'MagNum');
magnumval=getcellstruct(r,'MagNum',magnumind,1,1);

imn=1;
maggroupind=1;

while imn<length(magnumval)
    
    eqcount=1;
    exitwhile=0;
    
    while magnumval(imn)==magnumval(imn+eqcount) && ~exitwhile
        if imn+eqcount<length(magnumval)
            eqcount=eqcount+1;
        else
            exitwhile=1;
        end
    end
    
    if exitwhile
        eqcount=eqcount+1;
    end
    
    maggroups{maggroupind}=magnumind(imn:imn+eqcount-1); %#ok<*AGROW>
    imn=imn+eqcount;
    maggroupind=maggroupind+1;
    
end

return