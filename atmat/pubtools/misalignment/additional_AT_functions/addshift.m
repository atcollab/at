function addshift(ELEMINDEX, DX, DY)
%ADDSHIFT adds misalignment vectors T1, T2 for elements
%
% ADDSHIFT(ELEMINDEX, DX, DY) sets the entrance and exit misalignment vectors
%  of one element or a group of elements in the globally defined lattice THERING.
%  
%  DX, DY are displacements of the ELEMENT 
%  so the coordinate transformation on the particle at entrance is
%	X  ->  X-DX
%   Y  ->  Y-DY
%  The elements to be modified are given by ELEMINDEX  
%	Previous stored values are overwritten. 
%
% See also ADDXROT, ADDYROT, ADDSROT
 
global THERING
numelems = length(ELEMINDEX);

if (numelems ~= length(DX)) | (numelems ~= length(DY))
   error('ELEMINDEX, DX, and DY must have the same number of elements');
end

for i = 1:length(ELEMINDEX)
   V = zeros(1,6);
   V(1) = -DX(i);
   V(3) = -DY(i);
   curr_T1 = THERING{ELEMINDEX(i)}.T1;
   curr_T2 = THERING{ELEMINDEX(i)}.T2;
   THERING{ELEMINDEX(i)}.T1 = curr_T1 + V;
   THERING{ELEMINDEX(i)}.T2 = curr_T2 - V;
end
