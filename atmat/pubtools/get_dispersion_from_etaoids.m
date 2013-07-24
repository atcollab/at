function disp = get_dispersion_from_etaoids(ring,refpts)
%get_dispersion_from_etaoids computes dispersion functions (x,px,y,py) at refpts
%using the etaoids (E. Forest terminology) which are computed from the one turn map

for j=1:length(refpts)
    k=refpts(j);
    if k==1
        newring=ring;
    else
        newring=[ring(k:end);ring(1:k-1)];
    end
    m66=findm66(newring);
    A=amat(m66);
    [H1,H2,H3]=find_etaoids(A);
    H3(6,6)
    disp(:,j) = H3*[0 0 0 0 1 0]';
end