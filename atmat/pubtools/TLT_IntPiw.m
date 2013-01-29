function I=TLT_IntPiw(u,um,B1,B2) 
% integral in Piwinski Formula for the Lifetime

%sqrt(pi*(B1.^2-B2.^2))*um*
I=(...
    (...
    ((2+1./u).^2).*(((u./um)./(1+u))-1) +...
    1 - ...
    sqrt((1+u))./sqrt((u./um)) -...
    1./(2.*u).*(4+(1./u)).*log((u./um)./(1+u))...
    ).*exp(-B1.*u).*besseli(0,B2.*u).*sqrt(u)./sqrt((1+u))...
    );

% % % avoid Nan
% indNan=isnan(I);
% I(indNan(1):end)=0;

I(isnan(I))=0;
I(isinf(I))=0;%max(I(~isinf(I))); %

% infind=isinf(I);
% find(infind)
% I(infind)=I(find(infind)+ones(size(find(infind))));

return