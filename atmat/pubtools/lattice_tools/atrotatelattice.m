function newring = atrotatelattice(ring,index)
%atrotatelattice gives a new lattice starting with index as the first
%element.  We assume that the lattice is a collumn cell array

if index>length(ring);
error('index too big!');
end

if index==1
    newring=ring;
else
newring=[ring(index:end);ring(1:index-1)];
end