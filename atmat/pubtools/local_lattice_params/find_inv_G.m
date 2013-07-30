function [G1 G2 G3]=find_inv_G(m66)

[V,d]=eigs(m66);
v1=V(:,1);
v2=V(:,2);
v3=V(:,3);
v4=V(:,4);
v5=V(:,5);
v6=V(:,6);

%compute norms of each:
for j=1:6
n(j)=-(1/2)*i*V(:,j)'*jmat(3)*V(:,j);
end

%put positive modes before negative modes (swap columns if first n in pair
%is negative)
for j=1:3
    if n(2*j-1) <0 
        V(:,[2*j-1,2*j]) = V(:,[2*j,2*j-1]);
        n([2*j-1,2*j]) = n([2*j,2*j-1]);
    end
end
%now normalize each vector

for j=1:6
    Vn(:,j)=V(:,j)/sqrt(abs(n(j)));
end

%find the vecs that project most onto x,y,z, and reorder
%call these vectors Vxyz, ordered via (vx,vmx,vy,vmy,vz,vmz)
Vxyz(:,1)=(1/sqrt(2))*[i 1 0 0 0 0]';
Vxyz(:,2)=i*conj(Vxyz(:,1));
Vxyz(:,3)=(1/sqrt(2))*[0 0 i 1 0 0]';
Vxyz(:,4)=i*conj(Vxyz(:,3));
Vxyz(:,5)=(1/sqrt(2))*[0 0 0 0 i 1]';
Vxyz(:,6)=i*conj(Vxyz(:,5));

for j=1:6
    for k=1:6
        nn(j,k) = (1/2)*abs(sqrt(-i*Vn(:,j)'*jmat(3)*Vxyz(:,k)));
    end
end

%now use nn to order pairs
n1 = nn(1,[1,3,5]);
n3 = nn(3,[1,3,5]);
n5 = nn(5,[1,3,5]);

[vals1,sind1]=sort(n1,'descend');
[vals3,sind3]=sort(n3,'descend');
[vals5,sind5]=sort(n5,'descend');

orderp = [2*sind1(1)-1,2*sind3(1)-1,2*sind5(1)-1];
ordern = [2*sind1(1),2*sind3(1),2*sind5(1)];

V_ordered = Vn;
V_ordered(:,[1,3,5])=V_ordered(:,orderp); %reorder positive modes
V_ordered(:,[2,4,6])=V_ordered(:,ordern); %reorder negative modes


v1=V_ordered(:,1);
v2=V_ordered(:,3);
v3=V_ordered(:,5);


G1 = -(1/2)*jmat(3)*(conj(v1)*v1.'+v1*v1')*jmat(3);
G2 = -(1/2)*jmat(3)*(conj(v2)*v2.'+v2*v2')*jmat(3);
G3 = -(1/2)*jmat(3)*(conj(v3)*v3.'+v3*v3')*jmat(3);