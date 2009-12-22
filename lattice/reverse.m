function z = reverse(A)
%REVERSE reverses the order of elements in a one-dimensional MATLAB ARRAY
%
% ARRAY2 = REVERSE(ARRAY1)
%   ARRAY1 can be numeric, char (string), or cell array
L = length(A);
if iscell(A)
   z = cell(size(A));
   for i=1:L
      z{i} = A{L-i+1};
   end
elseif isnumeric(A)
   z = zeros(size(A));
   for i=1:L
      z(i) = A(L-i+1);
   end
elseif ischar(A)
    z = char(reverse(double(A)));
end
