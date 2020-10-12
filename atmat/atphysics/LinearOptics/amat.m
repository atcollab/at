function a=amat(transmat)
%find A matrix from one turn map matrix T such that:
%
%           [Rotx  0    0  ]
%inv(A)*T*A=[ 0   Rotz  0  ]
%           [ 0    0   Rots]
%
%Order it so that it is close to the order of x,y,z
%also ensure that positive modes are before negative so that
%one has proper symplecticity
%B. Nash July 18, 2013
%we order and normalize the vectors via
% v_j' jmat(3) v_k = i sgn(j) delta(j,k)

persistent Vxyz jm
if isempty(Vxyz)
    Vxyz={...
        1/sqrt(2)*[-1i -1;1 1i;],...
        1/sqrt(2)*[-1i -1 0 0;1 1i 0 0;0 0 -1i -1;0 0 1 1i;],...
        1/sqrt(2)*[-1i -1 0 0 0 0;1 1i 0 0 0 0;0 0 -1i -1 0 0;...
        0 0 1 1i 0 0; 0 0 0 0 -1i -1;0 0 0 0 1 1i;]...
        };
    jm={jmat(1),jmat(2),jmat(3)};
end

nv=size(transmat,1);
dms=nv/2;
idx=reshape(1:nv,2,dms);
select=idx(1,:);

[V,~]=eig(transmat);

%compute norms of each:
Vp=V'*jm{dms};
n=-0.5i*sum(Vp.'.*V);

%put positive modes before negative modes (swap columns if first n in pair
%is negative)
order=reshape(idx([1 1],:),1,nv) + (n<0);
V=V(:,order);
n=n(order);

%now normalize each vector
Vn=V./repmat(sqrt(abs(n)),nv,1);

%find the vecs that project most onto x,y,z, and reorder
%nn will have structure
% n1x n1y n1z
% n2x n2y n2z
% n3x n3y n3z
nn=0.5*abs(sqrt(-1i*Vn'*jm{dms}*Vxyz{dms}));

rows=select;
V_ordered=[];
for ixz=select
    [~,ind]=max(nn(rows,ixz));
    V_ordered=[V_ordered Vn(:,rows(ind))]; %#ok<AGROW>
    rows(ind)=[];
end

%build a matrix
a=reshape([real(V_ordered);imag(V_ordered)],nv,nv);

end


