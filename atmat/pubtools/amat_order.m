function a=amat_order(m66)

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
vx=(1/sqrt(2))*[1 i 0 0 0 0]';
vy=(1/sqrt(2))*[0 0 1 i 0 0]';
vz=(1/sqrt(2))*[0 0 0 0 1 i]';

n1x=abs(sqrt(-i*v1'*jmat(3)*vx));
n1y=abs(sqrt(-i*v1'*jmat(3)*vy));
n1z=abs(sqrt(-i*v1'*jmat(3)*vz));

n2x=abs(sqrt(-i*v2'*jmat(3)*vx));
n2y=abs(sqrt(-i*v2'*jmat(3)*vy));
n2z=abs(sqrt(-i*v2'*jmat(3)*vz));

n3x=abs(sqrt(-i*v3'*jmat(3)*vx));
n3y=abs(sqrt(-i*v3'*jmat(3)*vy));
n3z=abs(sqrt(-i*v3'*jmat(3)*vz));

[vals1,sind1]=sort([n1x,n1y,n1z],'descend')
[vals2,sind2]=sort([n2x,n2y,n2z],'descend')
[vals3,sind3]=sort([n3x,n3y,n3z],'descend')

Vecs=[v1,v2,v3]
v1x=Vecs(:,sind1(1))
v2y=Vecs(:,sind2(1))
v3z=Vecs(:,sind3(1))

a=[real(v1x) imag(v1x) real(v2y) imag(v2y) real(v3z) imag(v3z)];