function re=setXshift(r,pos,DX)
%
%
% See also SETSHIFT
 
numelems = length(pos);

if (numelems ~= length(DX))
   error('ELEMINDEX, DX must have the same number of elements');
end

for i = 1:length(pos)
   r{pos(i)}.T1(1) =  -DX(i);
   r{pos(i)}.T2(1) =   DX(i);
end

   re=r;


