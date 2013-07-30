function newring = atsetenergy(ring,E)
%atsetenergy(ring,Energy) sets the Energy field in all
%elements to the value given by E.  If no such field exists, it creates it.
newring=ring;
for j=1:length(ring)
    newring{j}.Energy = E;
end