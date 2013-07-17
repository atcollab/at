function [G1 G2 G3]=find_inv_G(m66)

[V,d]=eigs(m66);
v1=V(:,1);
v2=V(:,3);
v3=V(:,5);

%now normalize
n1=sqrt(-i*v1'*jmat(3)*v1);
n2=sqrt(-i*v2'*jmat(3)*v2);
n3=sqrt(-i*v3'*jmat(3)*v3);

v1=sqrt(2)*v1/n1;
v2=sqrt(2)*v2/n2;
v3=sqrt(2)*v3/n3;

G1 = -(1/2)*jmat(3)*(conj(v1)*v1.'+v1*v1')*jmat(3);
G2 = -(1/2)*jmat(3)*(conj(v2)*v2.'+v2*v2')*jmat(3);
G3 = -(1/2)*jmat(3)*(conj(v3)*v3.'+v3*v3')*jmat(3);