function [beam1,beam2,tune1,tune2] = beam44(A,B,C,gamma)
%BEAM44 computes the coupled beam matrices
%
%[BEAM1,BEAM2]=BEAM44(A,B,C,GAMMA)
% A,B,C,gamma: Coupling parameters, see [1]
% BEAM1,BEAM2: Eigen modes
%
%[BEAM1,BEAM2]=BEAM44(LINDATA)
% LINDATA: structure with fields A,B,C,gamma
%
%[BEAM1,BEAM2,TUNE1,TUNE1]=BEAM44(...)
% also returns the tunes
%
%[1] Sagan, Rubin, "Linear Analysis of Coupled Lattices"
%    Phys.Rev.Spec.Top. - Accelerators and Beams, vol2, 1999

if isstruct(A)
   [beam1,beam2,tune1,tune2]=beam44(A.A,A.B,A.C,A.gamma);
else
   S2=[0 1;-1 0];
   g=[gamma 0;0 gamma];
   v=[g C;-S2*C'*S2' g];
   [beama,tune1]=beam22(A);
   [beamb,tune2]=beam22(B);
   beam1=v*[beama zeros(2,2);zeros(2,4)]*v';
   beam2=v*[zeros(2,4);zeros(2,2) beamb]*v';
end
end

