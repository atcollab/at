function a=amat(m66)
%find A matrix from one turn map matrix
%Order it so that it is close to the order of x,y,z
%also ensure that positive modes are before negative so that
%one has proper symplecticity
%B. Nash July 18, 2013
%we order and normalize the vectors via
% v_j' jmat(3) v_k = i sgn(j) delta(j,k)

persistent Vxyz jm
if isempty(Vxyz)
    jm=jmat(3);
    Vxyz=1/sqrt(2)*[...
        -1i -1 0 0 0 0;...
        1 1i 0 0 0 0;...
        0 0 -1i -1 0 0;...
        0 0 1 1i 0 0;...
        0 0 0 0 -1i -1;...
        0 0 0 0 1 1i;...
        ];
end

[V,~]=eig(m66);

%compute norms of each:
Vp=V'*jm;
n=-0.5*1i*sum(Vp.'.*V);

%put positive modes before negative modes (swap columns if first n in pair
%is negative)
order=[1 1 3 3 5 5]+(n<0);
V=V(:,order);
n=n(order);

%now normalize each vector
Vn=V./repmat(sqrt(abs(n)),6,1);

%find the vecs that project most onto x,y,z, and reorder
nn=0.5*abs(sqrt(-1i*Vn'*jm*Vxyz));
[~,ind]=max(nn([1 3 5],[1 3 5]));

%reorder pairs
V_ordered=Vn(:,2*ind-1);

%build a matrix
a=reshape([real(V_ordered);imag(V_ordered)],6,6);

end


