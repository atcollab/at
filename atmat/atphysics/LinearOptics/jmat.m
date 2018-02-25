function mat=jmat(dim)
% JMAT Compute antisymmetric Matrix [O 1; -1 0]
% 
%  INPUTS
%    1. dim - 1,2,3 Dimension of the sqare matrix
%
%  OUPUTS
%    2. mat - Antisymmetric block Matrix [O 1; -1 0] 
% 
%  See also symplectify

S2 = [0 1; -1 0];

if(dim==1)
    mat=[0 1;-1 0];
elseif(dim==2)
    mat = blkdiag(S2,S2);
elseif(dim==3)
    mat = blkdiag(S2,S2,S2);
else
    Error('Dim is 1, 2 or 3')
end