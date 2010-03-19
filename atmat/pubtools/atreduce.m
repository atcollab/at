function [newring,seq,newrefs] = atreduce(oldring,oldrefs)
%ATREDUCE Remove useless elements from an AT structure
%NEWRING=ATREDUCE(OLDRING)
%
% Remove elements with PassMethod='IdentityPass' and merges adjacent drifts
%
%[NEWRING,KEPT]=ATREDUCE(OLDRING)
%
% Returns the index of kept elements
%
%[NEWRING,KEPT,NEWREFPTS]=ATREDUCE(OLDRING,OLDREFPTS)
%
% Returns in addition the updated list of reference points. Reference
% points are kept intact.
%

refs=false(1,length(oldring(:))+1);
seq=1:length(oldring(:));
if nargin >= 2, refs(oldrefs)=true; end
newring=oldring(:);
%					Remove useless elements
keep=true(1,length(newring));
keep(findcells(newring,'PassMethod','IdentityPass'))=false;
keep=keep | refs(1:end-1);
newring=newring(keep);
seq=seq(keep);
refs(~keep)=[];
%					Merge adjacent drifts
drift=false(1,length(newring));
drift(findcells(newring,'PassMethod','DriftPass'))=true;
drift(refs(1:end-1))=false;
lg=zeros(size(drift));
lg(drift)=getcellstruct(newring,'Length',find(drift));
newdrift=drift;
while true
   del=[newdrift(2:end) false] & newdrift & ~[false newdrift(1:end-1)];
   if ~any(del), break; end
   nex=[false del(1:end-1)];
   lg(nex)=lg(nex)+lg(del);
   lg(del)=0;
   newdrift=newdrift & ~del;
end
keep=~drift | newdrift;
newring=setcellstruct(newring,'Length',find(drift),lg(drift));
newring=newring(keep);
seq=seq(keep);
refs(~keep)=[];

if nargout >= 3, newrefs=find(refs); end
end
