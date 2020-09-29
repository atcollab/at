function [bx,bz] = Baxis(wig,s)
%Baxis   Compute the field in the axis of a generic wiggler
%
%The field of an horizontal wiggler is represented by a sum of harmonics,
%each being described by a 6x1 column in a matrix By such as:
%
%   Bz/Bmax = -B2 * cos(B5*kw*s + B6)
%
%The field of an vertical wiggler is represented by a sum of harmonics,
%each being described by a 6x1 column in a matrix Bx such as:
%
%   Bx/Bmax =  B2 * cos(B5*kw*s + B6)

kw=2*pi/wig.Lw;
Bmax=wig.Bmax;
kws=kw*s;
[bxh,bzh]=cellfun(@harmh,num2cell(wig.By,1),'UniformOutput',false);
[bxv,bzv]=cellfun(@harmv,num2cell(wig.Bx,1),'UniformOutput',false);
bx=sum(cat(3,bxh{:},bxv{:}),3);
bz=sum(cat(3,bzh{:},bzv{:}),3);

    function [bx,bz]=harm(pb)
        bz=-Bmax*pb(2)*cos(pb(5)*kws+pb(6));
        bx=zeros(size(kws));
    end
    function [bx,bz]=harmh(pb)
        [bx,bz]=harm(pb);
    end
    function [bx,bz]=harmv(pb)
        [bz,fz]=harm(pb);
        bx=-fz;
    end
end
