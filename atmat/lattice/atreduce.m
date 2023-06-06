function [newring,indx,newrefs] = atreduce(oldring,oldrefs)
%ATREDUCE Remove useless elements from an AT structure
%NEWRING=ATREDUCE(RING)
%
% Remove elements with PassMethod='IdentityPass' and merges adjacent
% similar elements
%
%NEWRING=ATREDUCE(RING,REFPTS)
%	When merging similar elements, keep REFPTS intact.
%
%[NEWRING,KEPT]=ATREDUCE(...)
%	Returns the index of kept elements so that NEWRING=OLDRING(KEPT)
%
%[NEWRING,KEPT,NEWREFPTS]=ATREDUCE(RING,REFPTS)
%	Returns in addition the updated list of reference points.
%

lg=numel(oldring);
refs=false(lg+1,1);
if nargin < 2
elseif islogical(oldrefs)
    refs(1:length(oldrefs))=oldrefs(:);
else
    refs(oldrefs)=true;
end
refend=refs(lg+1);
refs=refs(1:lg);
%					Remove useless elements
keep=~atgetcells(oldring(:),'PassMethod','IdentityPass') | refs;
newring=oldring(keep);
indx=keep;
%					Merge adjacent elements
tobecompared={'PassMethod','PolynomA','PolynomB','MaxOrder','RApertures','EApertures'};
bring=cellfun(@elstrip,newring,'UniformOutput',false);
bring(refs(indx))={NaN};
bends=atgetcells(newring,'BendingAngle');
keep=true(size(newring));
ba=zeros(size(newring));
inangle=zeros(size(newring));
outangle=zeros(size(newring));
invrad=zeros(size(newring));
ll=atgetfieldvalues(newring,'Length');
ba(bends)=atgetfieldvalues(newring(bends),'BendingAngle');
inangle(bends)=atgetfieldvalues(newring(bends),'EntranceAngle');
outangle(bends)=atgetfieldvalues(newring(bends),'ExitAngle');
invrad(bends)=ba(bends)./ll(bends);
islikenext=cellfun(@isequal,[bring(2:end);{NaN}],bring) ...
    & abs(invrad([2:end 1])-invrad)<5*eps(invrad) ...
    & abs(inangle([2:end 1])+outangle)<1.e-6;
islikeprev=circshift(islikenext,1);
arrayfun(@group,find(islikenext & ~islikeprev),find(~islikenext & islikeprev));
newring=newring(keep);
indx(indx)=keep;

if nargout >= 3
    if nargin >= 2 && islogical(oldrefs)
        if length(oldrefs)==lg
            newrefs=refs(indx);
        else
            newrefs=[refs(indx);refend];
        end
    else
        newrefs=find([refs(indx);refend]);
    end
end

    function el2=elstrip(el)
        if isfield(el,'Length') && el.Length>0
            el2=mvfield(struct(),el,tobecompared);
        else
            el2=NaN;
        end
    end

    function group(i1,i2)   % group consecutive elements
        newelem=newring{i1};
        newelem.Length=sum(ll(i1:i2));
        if bends(i1)
            newelem.BendingAngle=sum(ba(i1:i2));
        end
        newelem.FamName=gname(atgetfieldvalues(newring(i1:i2),'FamName'));
        newring{i1}=mvfield(newelem,newring{i2},exitfields());
        keep(i1+1:i2)=false;
    end

    function nm=gname(names)    % generate a group name as the longest common part
        lmin=min(cellfun(@length,names));
        nm=names{1};
        cellfun(@cmps,names);
        if lmin > 0, nm=nm(1:lmin); end
        
        function cmps(name)
            lmin=find([name(1:lmin) 'a'] ~= [nm(1:lmin) 'b'],1)-1;
            if (lmin>0) && (nm(lmin)=='_' || nm(lmin)=='-'), lmin=lmin-1; end
        end
    end
end
