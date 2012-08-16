function ferrs = atmakefielderrstruct(ring,class,nval,Bval,Aval,radius,rnd)
% MAKERNDFIELDERRS will create a field error data structure
%
% class is an element class, "Quadrupole", or "Sextupole"

if isequal(class,'Quadrupole')
        elemindex=findcells(ring,'Class','Quadrupole'); 
       Nval=2*ones(length(elemindex),1);   
elseif isequal(class,'Sextupole')
         elemindex=findcells(ring,'Class','Sextupole');
        Nval=3*ones(length(elemindex),1); 
else
    error('not a valid element class');
end

len=length(elemindex);
if(rnd==1)
for j=1:len
    Barray(j)=Bval*randn;
    Aarray(j)=Aval*randn;
end
else
    Barray=Bval*ones(len,1);
    Aarray=Aval*ones(len,1);
end

%e.g. ferr=struct('elemindex',[1 5 6],'nval',[4 4 4],'Nval',[2 2 2],'Bval',[1e-4 -1e-4 3e-4],'Aval',[0 0 0])

ferrs=struct('elemindex',elemindex,'Nval',Nval,'nval',nval*ones(len,1),'Bval',Barray,'Aval',Aarray,'radius',radius);