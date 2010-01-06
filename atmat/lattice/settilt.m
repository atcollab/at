function settilt(ELEMINDEX, PSI)
%SETTILT sets the entrance and exit misalignment matrixes
% of an element or a group of elements in THERING
% Previously stored values are overwritten.
%
% SETTILT(ELEMINDEX, PSI) 
% ELEMINDEX contains indexes of elements to be rotated
% PSI - angle(s) of rotation in RADIANS
%   POSITIVE PSI corresponds to a CORKSCREW (right) 
%   rotation of the ELEMENT.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The misalgnment matrixes are stored in fields R1 and R2
%   R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
%   R2 = R1'
%
% See also SETSHIFT MKSROLLMAT

global THERING

numelems = length(ELEMINDEX);

if numelems ~= length(PSI)
   error('ELEMINDEX and PSI must have the same number of elements');
end

C = cos(PSI);
S = sin(PSI);

for i = 1:length(ELEMINDEX)
   RM = diag([ C(i) C(i) C(i) C(i) 1  1 ]);
   RM(1,3) = S(i);
   RM(2,4) = S(i);
   RM(3,1) = -S(i);
   RM(4,2) = -S(i);
   THERING{ELEMINDEX(i)}.R1 = RM;
   THERING{ELEMINDEX(i)}.R2 = RM';
end

   