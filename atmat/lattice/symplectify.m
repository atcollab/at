function MS=symplectify(M)
%symplectify makes a matrix more symplectic
%follow Healy algorithm as described by McKay
%BNL-75461-2006-CP

J=jmat(3);
V=J*(eye(6)-M)*inv(eye(6)+M);
%V should be almost symmetric.  Replace with symmetrized version.

W=(V+V')/2;
%Now reconstruct M from W
MS=(eye(6)+J*W)*inv(eye(6)-J*W);