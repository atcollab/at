function [newring,seq,newrefs] = atreduce(oldring,oldrefs)
%ATREDUCE Remove useless elements from an AT structure
%NEWRING=ATREDUCE(OLDRING)
%
% Remove elements with PassMethod='IdentityPass' and merges adjacent
% similar elements
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
seq=1:lg;
%					Remove useless elements
keep=~atgetcells(oldring(:),'PassMethod','IdentityPass') | refs;
newring=oldring(keep);
seq=seq(keep);
refs=refs(keep);
%					Merge adjacent elements
bring=cellfun(@(el) elstrip(el,{'PassMethod','PolynomA','PolynomB','K','MaxOrder','NumIntSteps'}),...
    newring,'UniformOutput',false);
bends=atgetcells(newring,'BendingAngle');
keep=true(size(newring));
ba=zeros(size(newring));
invrad=zeros(size(newring));
ll=atgetfieldvalues(newring,'Length');
ba(bends)=atgetfieldvalues(newring(bends),'BendingAngle');
invrad(bends)=ba(bends)./ll(bends);
islikenext=cellfun(@isequal,[bring(2:end);{NaN}],bring) & abs(invrad([2:end 1])-invrad)<5*eps(invrad);
islikeprev=islikenext([end 1:end-1]);
arrayfun(@group,find(islikenext & ~islikeprev),find(~islikenext & islikeprev));
newring=newring(keep);
seq=seq(keep);
refs=refs(keep);

if nargout >= 3
    if nargin >= 2 && islogical(oldrefs)
        if length(oldrefs)==lg
            newrefs=refs;
        else
            newrefs=[refs;refend];
        end
    else
        newrefs=find([refs;refend]);
    end
end

    function el2=elstrip(el,fname)
        if isfield(el,'Length') && el.Length>0
            fnms=fieldnames(el);
            dd=~cellfun(@(ff) any(strcmp(ff,fname)),fnms);
            el2=rmfield(el,fnms(dd));
        else
            el2=NaN;
        end
    end

    function group(i1,i2)
        newelem=newring{i1};
        newelem.Length=sum(ll(i1:i2));
        if bends(i1)
            newelem.BendingAngle=sum(ba(i1:i2));
            newelem.ExitAngle=newring{i2}.ExitAngle;
        end
        newelem.FamName=gname(atgetfieldvalues(newring(i1:i2),'FamName'));
        newring{i1}=newelem;
        keep(i1+1:i2)=false;
    end

    function nm=gname(names)
        lmin=min(cellfun(@length,names));
        nm=names{1};
        cellfun(@cmps,names);
        if lmin > 0, nm=nm(1:lmin); end
        
        function cmps(name)
            lmin=find([name(1:lmin) 'a'] ~= [nm(1:lmin) 'b'],1)-1;
            if nm(lmin) == '_' || nm(lmin) == '-', lmin=lmin-1; end
        end
    end
end
