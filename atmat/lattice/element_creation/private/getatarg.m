function [rsrc,resval]=getatarg(rsrc,resval,resname)
%GETATARG extracts a resource from a resource list
%[RSRCLIST,RSRCVAL]=GETATARG(RSRCLIST,DEFAULTVAL,RSRCNAME)

ik=find(strcmp(resname,rsrc(1:2:end-1)),1,'last');
if ~isempty(ik)
    resval=rsrc{2*ik};
    rsrc(2*ik-1:2*ik)=[];
end
end

