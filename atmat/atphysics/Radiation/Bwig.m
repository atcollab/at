function [bx,bz,bs] = Bwig(wig,x,z,s)
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
        [bz,bb,bs]=harm(pb,kwz,-kwx,kws);
        bx=-bb;
    end
end
