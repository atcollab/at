function re=setANYshift(r,pos,coord,D)
%SETANYSHIFT Adds to the existing shift errors additional D
% 
%  See also SETSHIFT

numelems = length(pos);

if coord<1 | coord>6
    error('coord should be an index vector between 1 and 6')
end

if (numelems ~= length(D))
    error('pos, D must have the same number of elements');
end

additive=0;% errors are overwritten if this flag is 0.

for i = 1:length(pos)
    if abs(D(i))<1e-12, D(i)=0; end;
    
    if additive
        r{pos(i)}.T1(coord,1) = r{pos(i)}.T1(coord,1) -D(i);
        r{pos(i)}.T2(coord,1) = r{pos(i)}.T1(coord,1) +D(i);
    else % no T1 and T2 defined
        r{pos(i)}.T1(coord,1) =  -D(i);
        r{pos(i)}.T2(coord,1) =   D(i);
    end
end

re=r;
