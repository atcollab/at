function I=TLT_IntPiw_k(k,km,B1,B2) 
% integral in Piwinski Formula for the Lifetime with u=tan^2(k)

t=tan(k).^2;
tm=tan(km).^2;

%   in case the Bessel function has too large value (more than 10^251) it 
%   is substituted by its exponential approximation:
%   I_0(x)~exp(x)/sqrt(2*pi*x)

if B2*t<500
    I=( ( (2.*t+1) .^2) .*( (t./tm)./(1+t) -1) ./t...
        +t...
        -sqrt(t.*tm.*(1+t))...
        -(2+1./(2.*t)).*log((t./tm)./(1+t)) ).*exp(-B1.*t).*besseli(0,B2.*t).*sqrt((1+t));
else
    I=( ( (2.*t+1) .^2) .*( (t./tm)./(1+t) -1) ./t...
        +t...
        -sqrt(t.*tm.*(1+t))...
        -(2+1./(2.*t)).*log((t./tm)./(1+t)) ).*exp(B2.*t-B1.*t)./sqrt(2*pi*B2.*t).*sqrt((1+t));
end

return