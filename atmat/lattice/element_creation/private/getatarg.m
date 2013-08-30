function [rsrc,resval]=getatarg(rsrc,resval,resname)
%GETATARG extracts a resource from a resource list
%[RSRCLIST,RSRCVAL]=GETATARG(RSRCLIST,DEFAULTVAL,RSRCNAME)

ok=strcmp(resname,rsrc(1:2:end-1));
if any(ok)
    resval=rsrc{2*find(ok,1,'last')};
    rsrc(reshape([ok;ok],1,[]))=[];
end
end

