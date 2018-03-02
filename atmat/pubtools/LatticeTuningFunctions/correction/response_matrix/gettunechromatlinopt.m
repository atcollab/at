function tc=gettunechromatlinopt(r,inCOD)
%GETTUNECHROMATLINOPT Gets tunes and chromaticities from atlinopt

[~,t,c]=atlinopt(r,0,1,inCOD);
tc=[t';c'];
end