function tc=gettunechromatlinopt(r,inCOD)
[~,t,c]=atlinopt(r,0,1,inCOD);
tc=[t';c'];
end