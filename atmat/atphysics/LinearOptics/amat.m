function [a,lambda]=amat(transmat)
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

nv=size(transmat,1);
dms=nv/2;
slices=num2cell(reshape(1:nv,2,dms),1);
Vxyz=makeVxyz(slices);
S=kmat(dms);
select=1:2:nv;

[V,lambda]=eig(transmat);
lambda=diag(lambda).';

%compute norms of each:
Vp=V'*S;
n=-0.5i*sum(Vp.'.*V);
if any(abs(n) < 1.0E-12)
    error('AT:Unstable','Unstable machine');
end

%put positive modes before negative modes (swap columns if first n in pair
%is negative)
order=reshape([select;select],1,nv) + (n<0);
V=V(:,order);
n=n(order);
lambda=lambda(order);

%now normalize each vector
%Vn=V./repmat(sqrt(abs(n)),nv,1);
Vn=V./sqrt(abs(n));

%find the vecs that project most onto x,y,z, and reorder
%nn will have structure
% n1x n1y n1z
% n2x n2y n2z
% n3x n3y n3z
nn=0.5*abs(sqrt(-1i*Vn'*S*Vxyz));

rows=select;
order=[];
for ixz=select
    [~,ind]=max(nn(rows,ixz));
    order=[order rows(ind)]; %#ok<AGROW>
    rows(ind)=[];
end
V_ordered=Vn(:,order);
lambda=lambda(order);

%build the a-matrix
a=reshape([real(V_ordered);imag(V_ordered)],nv,nv);

    function Vxyz=makeVxyz(slices)
        Vxyz=zeros(nv,nv);
        cellfun(@s22, slices)
        
        function s22(slc)
            Vxyz(slc,slc)=[-1i -1;1 1i];
        end
    end

    function mat=kmat(dim)
        % Modified version of jmat to deal with the swap of the
        % longitudinal coordinates
        S2 = [0 1; -1 0];
        
        if(dim==1)
            mat=S2;
        elseif(dim==2)
            mat = blkdiag(S2,S2);
        elseif(dim==3)
            mat = blkdiag(S2,S2,S2');
        else
            Error('Dim is 1, 2 or 3')
        end
    end

end


