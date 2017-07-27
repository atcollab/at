function RM = mkSRotationMatrix(psi);
%MKSROTATIONMATRIX(PSI) coordinate transformation matrix 
% that describes s-rotation of the ELEMENT 
% 
% |   cos(psi)     0       sin(psi)     0         0       0     |
% |      0     cos(psi)        0      sin(psi)    0       0     |
% |  -sin(psi)     0       cos(psi)     0         0       0     |
% |      0     -sin(psi)       0      cos(psi)    0       0     |
% |      0         0           0        0         1       0     |
% |      0         0           0        0         0       1     |
%
% Note: AT defines 'positive' s-rotation as clockwise,  
%       looking in the dirction of the beamm
%

C = cos(psi);
S = sin(psi);

RM = diag([ C C C C 1  1 ]);
RM(1,3) =  S;
RM(2,4) =  S;
RM(3,1) = -S;
RM(4,2) = -S;