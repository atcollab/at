function a=amat(m66)

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

%find the vecs that project most onto x,y,z, and reorder
vx=(1/sqrt(2))*[1 i 0 0 0 0 0]';
vy=(1/sqrt(2))*[0 0 1 i 0 0 0]';
vz=(1/sqrt(2))*[0 0 0 0 0 1 i]';

n1x=sqrt(-i*v1'*jmat(3)*vx);
n1y=sqrt(-i*v1'*jmat(3)*vy);
n1z=sqrt(-i*v1'*jmat(3)*vz);

n2x=sqrt(-i*v2'*jmat(3)*vx);
n2y=sqrt(-i*v2'*jmat(3)*vy);
n2z=sqrt(-i*v2'*jmat(3)*vz);

n3x=sqrt(-i*v3'*jmat(3)*vx);
n3y=sqrt(-i*v3'*jmat(3)*vy);
n3z=sqrt(-i*v3'*jmat(3)*vz);



a=[real(v1x) imag(v1x) real(v2y) imag(v2y) real(v3z) imag(v3z)];