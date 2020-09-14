function [bx,bz,bs] = Bwig(wig,x,z,s)
%Bwig   Compute the field of a geeric wiggler
%
%The field of an horizontal wiggler is represented by a sum of harmonics,
%each being described by a 6x1 column in a matrix By such as:
%
%   Bx/Bmax =  B2 * B3/B4 * sin(B3*kw*x) sinh(B4*kw*z) cos(B5*kw*s + B6)
%   Bz/Bmax = -B2 * cos(B3*kw*x) cosh(B4*kw*z) cos(B5*kw*s + B6)
%   Bs/Bmax =  B2 * B5/B4 * cos(B3*kw*x) sinh(B4*kw*z) sin(B5*kw*s + B6)
%
%The field of an vertical wiggler is represented by a sum of harmonics,
%each being described by a 6x1 column in a matrix Bx such as:
%
%   Bx/Bmax =  B2 * cos(B4*kw*z) cosh(B3*kw*x) cos(B5*kw*s + B6)
%   Bz/Bmax = -B2 * B4/B3 * sin(B4*kw*z) sinh(B3*kw*x) cos(B5*kw*s + B6)
%   Bs/Bmax = -B2 * B5/B3 * cos(B3*kw*z) sinh(B3*kw*x) sin(B5*kw*s + B6)

kw=2*pi/wig.Lw;
Bmax=wig.Bmax;
kwx=kw*x;
kwz=kw*z;
kws=kw*s;
[bxh,bzh,bsh]=cellfun(@harmh,num2cell(wig.By,1),'UniformOutput',false);
[bxv,bzv,bsv]=cellfun(@harmv,num2cell(wig.Bx,1),'UniformOutput',false);
bx=squeeze(sum(cat(4,bxh{:},bxv{:}),4));
bz=squeeze(sum(cat(4,bzh{:},bzv{:}),4));
bs=squeeze(sum(cat(4,bsh{:},bsv{:}),4));

    function [bx,bz,bs]=harm(pb,kwx,kwz,kws)
        coef=Bmax*pb(2);
        kx=pb(3)*kwx;
        kz=pb(4)*kwz;
        ks=pb(5)*kws+pb(6);
        sx=sin(kx);
        cx=cos(kx);
        shz=sinh(kz);
        chz=cosh(kz);
        cs=cos(ks);
        ss=sin(ks);
        [a,b,c]=ndgrid(sx,shz,cs);
        bx=coef*pb(3)/pb(4)*(a.*b.*c);
        [a,b,c]=ndgrid(cx,chz,cs);
        bz=-coef*(a.*b.*c);
        [a,b,c]=ndgrid(cx,shz,ss);
        bs=coef*pb(5)/pb(4)*(a.*b.*c);
    end
    function [bx,bz,bs]=harmh(pb)
        [bx,bz,bs]=harm(pb,kwx,kwz,kws);
    end
    function [bx,bz,bs]=harmv(pb)
        [fx,fz,fs]=harm(pb([1 2 4 3 5 6]),kwz,-kwx,kws);
        bx=-permute(fz,[2 1 3]);
        bz=permute(fx,[2 1 3]);
        bs=permute(fs,[2 1 3]);
    end
end
