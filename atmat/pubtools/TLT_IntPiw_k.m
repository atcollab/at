function I=TLT_IntPiw_k(k,km,B1,B2) 
% integral in Piwinski Formula for the Lifetime with u=tan^2(k)

t=tan(k).^2;
tm=tan(km).^2;

I=( ( (2.*t+1) .^2) .*( (t./tm)./(1+t) -1) ./t...
    +t...
    -sqrt(t.*tm.*(1+t))...
    -(2+1./(2.*t)).*log((t./tm)./(1+t)) ).*exp(-B1.*t).*besseli(0,B2.*t).*sqrt((1+t));

I(isnan(I))=0;
I(isinf(I))=0;

return